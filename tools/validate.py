#!/usr/bin/env python
"""Validate the pipeline on simulated nanopore data with known editing rates.

Simulates an unedited control plus edited samples at several editing levels
for two error regimes (R10.4.1-like and R9.4.1-like), runs the full CLI on
each regime, and compares measured vs true % modified / % frameshift.

  python tools/validate.py --out validation_run --reads 5000
Writes <out>/validation.tsv and <out>/validation.png
"""
import argparse
import csv
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)


def truth(prefix):
    rows = [r for r in csv.DictReader(open(prefix + "_truth.tsv"), delimiter="\t") if r["edited"] in "01"]
    n = len(rows)
    return (100 * sum(r["edited"] == "1" for r in rows) / n,
            100 * sum(r["category"] == "frameshift" for r in rows) / n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--locus", default=os.path.join(ROOT, "config", "PROCR_g1.yaml"))
    ap.add_argument("--out", default="validation_run")
    ap.add_argument("--reads", type=int, default=5000)
    ap.add_argument("--levels", default="0.1,0.25,0.5,0.75,0.9")
    ap.add_argument("--engine", default="both")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    regimes = {"R10-like (scale 1)": 1.0, "R9-like (scale 4)": 4.0}
    levels = [float(x) for x in a.levels.split(",")]
    out_rows = []
    for ri, (reg, scale) in enumerate(regimes.items()):
        d = os.path.join(a.out, f"sim_scale{scale:g}")
        os.makedirs(d, exist_ok=True)
        pref = {}
        for i, e in enumerate([0.0] + levels):
            p = os.path.join(d, f"edit{int(round(e * 100)):02d}")
            subprocess.run([sys.executable, os.path.join(HERE, "simulate_reads.py"), "--locus", a.locus,
                            "--platform", "ont", "--edit", str(e), "--pairs", str(a.reads),
                            "--error-scale", str(scale), "--seed", str(100 * ri + i + 1), "--out", p],
                           check=True, stdout=subprocess.DEVNULL)
            pref[e] = p
        res = os.path.join(d, "results")
        cmd = [sys.executable, os.path.join(ROOT, "run_crispr_indel.py"), "--locus", a.locus,
               "--control", pref[0.0] + "_R1.fastq.gz", "-o", res, "--engine", a.engine] + \
              [pref[e] + "_R1.fastq.gz" for e in levels]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        met = list(csv.DictReader(open(os.path.join(res, "all_metrics.tsv")), delimiter="\t"))
        for e in [0.0] + levels:
            sname = os.path.basename(pref[e]) + "_R1"
            tm, tf = truth(pref[e])
            for r in met:
                if r["sample"] == sname and r["window"] == "1":
                    out_rows.append(dict(regime=reg, true_pct_modified=tm, true_pct_frameshift=tf,
                                         engine=r["engine"], pct_modified=float(r["pct_modified"]),
                                         pct_modified_corrected=float(r["pct_modified_corrected"]),
                                         pct_frameshift=float(r["pct_frameshift"]),
                                         pct_frameshift_corrected=float(r["pct_frameshift_corrected"]),
                                         informative_reads=int(r["informative_reads"])))
    cols = list(out_rows[0].keys())
    with open(os.path.join(a.out, "validation.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out_rows:
            fh.write("\t".join(f"{r[c]:.3f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    print(f"wrote {a.out}/validation.tsv ({len(out_rows)} rows)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), sharex=True, sharey=True)
    for ax, reg in zip(axes, regimes):
        for eng, mk, col in (("CRISPResso2", "o", "#b2182b"), ("builtin", "s", "#2166ac")):
            rr = [r for r in out_rows if r["regime"] == reg and r["engine"] == eng]
            if not rr:
                continue
            ax.plot([r["true_pct_modified"] for r in rr], [r["pct_modified_corrected"] for r in rr], mk,
                    color=col, mfc=col, label=f"{eng}, control-corrected % modified", ms=6)
            ax.plot([r["true_pct_frameshift"] for r in rr], [r["pct_frameshift_corrected"] for r in rr], mk,
                    color=col, mfc="none", label=f"{eng}, corrected % frameshift", ms=6)
        ax.plot([0, 100], [0, 100], color="0.6", lw=1, ls="--")
        ax.set_title(reg, fontsize=10)
        ax.set_xlabel("true %")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("measured %")
    axes[0].legend(frameon=False, fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "validation.png"), dpi=200)
    print(f"wrote {a.out}/validation.png")


if __name__ == "__main__":
    main()
