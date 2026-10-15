"""crispr-indel: CRISPR nuclease editing efficiency from amplicon FASTQ.

Typical use (Oxford Nanopore, one FASTQ per sample, unedited control):

  crispr-indel --locus config/PROCR_g1.yaml --control WT.fastq.gz \\
               KO_pool1.fastq.gz KO_pool2.fastq.gz -o results/

Positional FASTQs may also be directories (all *.fastq / *.fastq.gz / *.fq.gz
inside are used).  For paired Illumina data give a samples TSV
(columns: sample, fastq_r1, fastq_r2, control[yes/no]) with --samples.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
import time
from collections import Counter
from pathlib import Path

from . import __version__
from . import crispresso as CR
from . import quick as Q
from . import report as R
from .locus import load_locus
from .preflight import preflight, sniff
from .prep import prep_fastq

FASTQ_RE = re.compile(r"\.(fastq|fq)(\.gz)?$")


def sample_name(path: str) -> str:
    b = os.path.basename(path)
    b = FASTQ_RE.sub("", b)
    return re.sub(r"[^A-Za-z0-9_.-]", "_", b)


def expand(paths):
    out = []
    for p in paths:
        if os.path.isdir(p):
            out += sorted(str(x) for x in Path(p).iterdir() if FASTQ_RE.search(x.name))
        else:
            out.append(p)
    return out


def detect_platform(r1, r2, locus, n=2000) -> str:
    if r2:
        return "illumina"
    s = sniff(r1, locus, n)
    # Nanopore/PacBio reads span the whole amplicon (plus flanks); Illumina SE reads are <= 300 nt
    if s.median_len > 320 or s.max_len > 600 or s.median_len >= 0.9 * locus.length:
        return "ont"
    return "illumina"


def top_alleles(allele_rows, k=15):
    agg = Counter()
    cat = {}
    tot = 0
    for lab, c, net, n, _ in allele_rows:
        agg[lab] += n
        cat[lab] = (c, net)
        tot += n
    return [(lab, cat[lab][0], cat[lab][1], n, 100 * n / tot) for lab, n in agg.most_common(k)]


def quick_allele_rows(calls):
    rows = []
    for _, st, c in calls:
        if st == "ok":
            rows.append((c.alleles or "WT", c.category, c.net, 1, ""))
    return rows


def analyze(sample, r1, r2, locus, a, outdir, log):
    sd = outdir / sample
    sd.mkdir(parents=True, exist_ok=True)
    platform = a.platform if a.platform != "auto" else detect_platform(r1, r2, locus)
    info = dict(sample=sample, fastq_r1=r1, fastq_r2=r2 or "", platform=platform)
    log(f"\n== {sample} ({platform}) ==")
    rows = []
    allele_rows_head = None
    if platform == "ont":
        st = sniff(r1, locus, 5000)
        log("  " + st.line())
        prepped = str(sd / "amplicon_reads.fastq.gz")
        pc = prep_fastq(r1, prepped, locus, min_q=a.min_q, max_reads=a.max_reads or None,
                        max_insertion=a.max_insertion)
        tot = sum(v for k, v in pc.items() if not k.startswith("strand"))
        info.update(reads_input=tot, reads_amplicon=pc.get("ok", 0), prep=dict(pc))
        log(f"  primer-anchored extraction: {dict(pc)}")
        if tot and pc.get("ok", 0) / tot < 0.5:
            log("  WARNING: <50% of reads carry both primers - check primers/sample identity/read quality.")
        fq1, fq2, oriented = prepped, None, True
    else:
        pf = preflight(locus, r1, r2, a.window)
        for m in pf.messages:
            log("  " + m)
        if not pf.ok and a.mode != "single":
            raise SystemExit(f"[{sample}] pre-flight failed; see messages above")
        info.update(mode=pf.mode)
        fq1, fq2, oriented = r1, (r2 if pf.mode == "merged" else None), False

    windows = sorted(set([a.window, 1, 5, 10, 0]), key=lambda w: (w != a.window, w if w else 999))
    # ---- engine 1: CRISPResso2
    cr_ok = a.engine in ("both", "crispresso") and CR.crispresso_available()
    if a.engine in ("both", "crispresso") and not cr_ok:
        log("  CRISPResso2 not found on PATH - using the built-in engine only (install: see README).")
    if cr_ok:
        t = time.time()
        cmd = CR.build_command(locus, fq1, fq2, sample, str(sd), a.window, platform, a.threads,
                               min_avg_q=0 if platform == "ont" else a.illumina_min_q,
                               exclude_bp=a.exclude_bp, count_substitutions=a.count_substitutions)
        (sd / "crispresso_command.txt").write_text(" ".join(cmd) + "\n")
        out = CR.run(cmd, sd / "crispresso.log")
        q = CR.read_quant(out)
        res, arows = CR.rescore(out, locus, windows, a.flank, a.exclude_bp, a.count_substitutions)
        allele_rows_head = arows
        for w in windows:
            m = R.metrics(res[w])
            rows.append(dict(sample=sample, engine="CRISPResso2", window=w, aligned_reads=res[w]["aligned"],
                             crispresso_native_pct_modified=float(q.get("Modified%", "nan")) if w == a.window else float("nan"),
                             **m))
        log(f"  CRISPResso2 {time.time() - t:.0f}s: native Modified% (window +/-{a.window}) = "
            f"{float(q.get('Modified%', 'nan')):.2f}; reads aligned {q.get('Reads_aligned')}")
    # ---- engine 2: built-in aligner (always runs; cross-check or fallback)
    if a.engine in ("both", "quick") or not cr_ok:
        t = time.time()
        if fq2 and platform == "illumina":
            qr = Q.run_paired(fq1, fq2, locus, a.window, a.flank, (a.exclude_bp, a.exclude_bp), a.max_reads or None)
        else:
            qr = Q.run_single(fq1, locus, a.window, a.flank, oriented, (a.exclude_bp, a.exclude_bp),
                              a.max_reads or None, platform=platform)
        m = R.metrics(qr.counts)
        rows.append(dict(sample=sample, engine="builtin", window=a.window,
                         aligned_reads=qr.counts["reads"] - qr.counts["unaligned"],
                         crispresso_native_pct_modified=float("nan"), **m))
        with open(sd / "builtin_read_calls.tsv", "w") as fh:
            fh.write("read\tstatus\tcategory\tnet\talleles\n")
            for n, st_, c in qr.calls:
                fh.write(f"{n}\t{st_}\t{c.category if c else ''}\t{c.net if c else ''}\t{c.alleles if c else ''}\n")
        if allele_rows_head is None:
            allele_rows_head = quick_allele_rows(qr.calls)
        log(f"  built-in engine {time.time() - t:.0f}s: {m['pct_modified']:.2f}% modified "
            f"({m['informative_reads']} informative reads)")
    tops = top_alleles(allele_rows_head or [])
    with open(sd / "top_alleles.tsv", "w") as fh:
        fh.write("allele\tcategory\tnet_indel\treads\tpct\n")
        for r in tops:
            fh.write("\t".join(map(str, r)) + "\n")
    return info, rows, allele_rows_head or [], tops


def main(argv=None):
    ap = argparse.ArgumentParser(prog="crispr-indel", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fastq", nargs="*", help="edited-sample FASTQ(s) or directories")
    ap.add_argument("--control", action="append", default=[],
                    help="unedited control FASTQ (strongly recommended for ONT); may repeat")
    ap.add_argument("--samples", help="TSV: sample, fastq_r1, [fastq_r2], [control yes/no]")
    ap.add_argument("--locus", help="locus YAML (amplicon, guide, primers, coding_seq)")
    ap.add_argument("--amplicon"), ap.add_argument("--guide"), ap.add_argument("--coding-seq")
    ap.add_argument("--forward-primer"), ap.add_argument("--reverse-primer")
    ap.add_argument("-o", "--out", default="crispr_indel_results")
    ap.add_argument("--platform", choices=["auto", "ont", "illumina"], default="auto")
    ap.add_argument("--engine", choices=["both", "crispresso", "quick"], default="both")
    ap.add_argument("-w", "--window", type=int, default=1,
                    help="quantification half-window around the cut (bp); 1 = CRISPResso2 default")
    ap.add_argument("--count-substitutions", action="store_true",
                    help="also count substitutions inside the window as edits (CRISPResso2 engine only). Recovers "
                         "some ONT reads where an error masks an indel, but raises background and counts real SNPs")
    ap.add_argument("--min-q", type=float, default=10.0,
                    help="ONT: minimum mean Phred quality of the primer-to-primer read (default 10)")
    ap.add_argument("--illumina-min-q", type=int, default=0, help="Illumina: CRISPResso2 --min_average_read_quality (CRISPResso2 default 0)")
    ap.add_argument("--max-insertion", type=int, default=150, help="ONT: discard products longer than amplicon+N")
    ap.add_argument("--exclude-bp", type=int, default=20, help="ignore indels within N bp of amplicon ends (primers)")
    ap.add_argument("--flank", type=int, default=15, help="read must extend N bp beyond the window to be informative")
    ap.add_argument("--mode", choices=["auto", "single"], default="auto", help="Illumina: force single-read use")
    ap.add_argument("--max-reads", type=int, default=0, help="analyse only the first N reads per sample (0 = all)")
    ap.add_argument("-p", "--threads", type=int, default=1, help="CRISPResso2 processes")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    a = ap.parse_args(argv)

    locus = load_locus(a.locus, a.amplicon, a.guide, a.coding_seq, a.forward_primer, a.reverse_primer)
    outdir = Path(a.out)
    outdir.mkdir(parents=True, exist_ok=True)
    logf = open(outdir / "run.log", "w")

    def log(m):
        print(m, flush=True)
        logf.write(m + "\n")

    log(f"crispr-indel {__version__}")
    log(locus.summary())
    for w in locus.warnings:
        log("WARNING: " + w)

    samples = []  # (name, r1, r2, is_control)
    if a.samples:
        for r in csv.DictReader(open(a.samples), delimiter="\t"):
            samples.append((r["sample"], r["fastq_r1"], r.get("fastq_r2") or None,
                            str(r.get("control", "")).lower() in ("yes", "y", "true", "1")))
    for p in expand(a.control):
        samples.append((sample_name(p), p, None, True))
    for p in expand(a.fastq):
        samples.append((sample_name(p), p, None, False))
    if not samples:
        ap.error("no input FASTQ given")
    if not any(s[3] for s in samples):
        log("NOTE: no --control given. Editing % will not be background-corrected; with nanopore reads the "
            "background at the cut site is typically 0.5-3%.")

    all_rows, infos, sizes = [], [], []
    alle_tables = {}
    for name, r1, r2, is_ctrl in samples:
        info, rows, arows, tops = analyze(name, r1, r2, locus, a, outdir, log)
        info["control"] = is_ctrl
        for r in rows:
            r["role"] = "control" if is_ctrl else "edited"
        infos.append(info)
        all_rows += rows
        alle_tables[name] = tops[:10]
        sizes.append((name, arows))

    # control correction, per engine/window: mean of control samples
    ctrl = {}
    for r in all_rows:
        if r["role"] == "control":
            ctrl.setdefault((r["engine"], r["window"]), []).append(r)
    for r in all_rows:
        cs = ctrl.get((r["engine"], r["window"]))
        if cs:
            cm = sum(c["pct_modified"] for c in cs) / len(cs)
            cf = sum(c["pct_frameshift"] for c in cs) / len(cs)
            r["control_pct_modified"] = cm
            r["pct_modified_corrected"] = R.correct(r["pct_modified"], cm)
            # mixture model: observed = E*true + (1-E)*control  ->  E*f = T_fs - (1-E)*C_fs
            e = r["pct_modified_corrected"] / 100
            r["pct_frameshift_corrected"] = max(0.0, r["pct_frameshift"] - (1 - e) * cf) \
                if not math.isnan(e) else float("nan")
        else:
            r["control_pct_modified"] = r["pct_modified_corrected"] = r["pct_frameshift_corrected"] = float("nan")

    cols = ["sample", "role", "engine", "window", "aligned_reads", "informative_reads", "pct_modified",
            "pct_modified_ci95_lo", "pct_modified_ci95_hi", "pct_frameshift", "pct_inframe", "pct_exon_boundary",
            "pct_noncoding_indel", "control_pct_modified", "pct_modified_corrected", "pct_frameshift_corrected",
            "crispresso_native_pct_modified"]
    with open(outdir / "all_metrics.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in all_rows:
            fh.write("\t".join(f"{r.get(c):.3f}" if isinstance(r.get(c), float) else str(r.get(c, ""))
                               for c in cols) + "\n")
    # headline: CRISPResso2 at the chosen window if available, else built-in
    head = []
    for name, *_ in samples:
        cand = [r for r in all_rows if r["sample"] == name and r["window"] == a.window]
        pick = next((r for r in cand if r["engine"] == "CRISPResso2"), cand[0] if cand else None)
        other = next((r for r in cand if r is not pick), None)
        if pick:
            h = dict(pick)
            h["builtin_pct_modified"] = other["pct_modified"] if other else float("nan")
            if other and abs(other["pct_modified"] - pick["pct_modified"]) > 2 and not a.count_substitutions:
                log(f"WARNING [{name}]: engines disagree by {abs(other['pct_modified'] - pick['pct_modified']):.1f} "
                    "percentage points - inspect top_alleles.tsv / CRISPResso2 report.")
            head.append(h)
    hcols = ["sample", "role", "engine", "window", "informative_reads", "pct_modified", "pct_modified_ci95_lo",
             "pct_modified_ci95_hi", "pct_frameshift", "pct_inframe", "pct_exon_boundary", "control_pct_modified",
             "pct_modified_corrected", "pct_frameshift_corrected", "builtin_pct_modified"]
    with open(outdir / "editing_summary.tsv", "w") as fh:
        fh.write("\t".join(hcols) + "\n")
        for r in head:
            fh.write("\t".join(f"{r.get(c):.2f}" if isinstance(r.get(c), float) else str(r.get(c, ""))
                               for c in hcols) + "\n")
    (outdir / "samples.json").write_text(json.dumps(infos, indent=1, default=str))

    figs = [("Outcome composition (headline engine/window); error bars = Wilson 95% CI of % modified",
             R.plot_outcomes(head, str(outdir / "outcomes.png")))]
    for name, arows in sizes:
        if arows:
            figs.append((f"{name}: indel sizes", R.plot_sizes(name, arows, str(outdir / f"{name}_indel_sizes.png"))))
    notes = [
        "pct_modified: reads with an insertion/deletion overlapping the cut-site window (substitutions ignored) / "
        "informative reads (reads that span the window).",
        "pct_frameshift: net indel length within the coding sequence not divisible by 3 - the best read-level "
        "proxy for loss-of-function alleles; in-frame indels may or may not disrupt function.",
        "pct_*_corrected: background-corrected with the unedited control, (T - C)/(100 - C). For nanopore data "
        "this removes sequencing-error indels at the cut site and is the number to report.",
        "window: half-size in bp around the cut (1 = CRISPResso2 default, 0 = whole amplicon minus primers).",
        "Large deletions that remove a primer site and alleles from CRISPR-induced large rearrangements are "
        "invisible to any amplicon assay.",
    ]
    pre = [locus.summary()] + [f"{i['sample']}: platform={i['platform']} "
                               + (f"reads={i.get('reads_input')} amplicon_reads={i.get('reads_amplicon')}"
                                  if i["platform"] == "ont" else f"mode={i.get('mode')}") for i in infos]
    R.html_report(str(outdir / "report.html"), f"CRISPR indel report - {locus.name}", pre, head, hcols, figs,
                  alle_tables, notes)
    log(f"\nWrote {outdir}/editing_summary.tsv, all_metrics.tsv, report.html")
    for r in head:
        corr = r.get("pct_modified_corrected")
        log(f"  {r['sample']:<30} {r['role']:<8} modified {r['pct_modified']:6.2f}% "
            f"[{r['pct_modified_ci95_lo']:.1f}-{r['pct_modified_ci95_hi']:.1f}]  frameshift {r['pct_frameshift']:6.2f}%"
            + (f"  corrected {corr:6.2f}%" if isinstance(corr, float) and not math.isnan(corr) else ""))
    logf.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
