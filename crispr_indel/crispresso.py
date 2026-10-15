"""Run CRISPResso2 and re-score its allele table.

CRISPResso2 (Clement et al., Nat Biotechnol 2019) is the reference engine.  We
call it with settings appropriate for SpCas9 indel quantification and then
re-score its full Alleles_frequency_table with crispr_indel.indels so that:
  * several quantification windows (+/-1, 5, 10 bp, whole amplicon) are
    reported from one run,
  * reads that do not span the window (terminal gaps) are excluded rather than
    counted as unmodified,
  * frameshift / in-frame calls use the exon-3 coding sequence and are defined
    identically for both engines.
CRISPResso2's own headline numbers are reported alongside.
"""
from __future__ import annotations

import csv
import io
import shutil
import subprocess
import zipfile
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional

from .indels import classify, parse_alignment
from .locus import Locus

# CRISPResso2 Needleman-Wunsch presets.  'ont' lowers the gap-open penalty so a
# nanopore sequencing indel near the cut does not cancel a true indel into
# mismatches (see docs/METHOD.md and the validation in tests/).
NW_PRESETS = {
    "illumina": dict(gap_open=-20, gap_extend=-2, gap_incentive=1),
    "ont": dict(gap_open=-10, gap_extend=-2, gap_incentive=1),
}


def crispresso_available() -> bool:
    return shutil.which("CRISPResso") is not None


def build_command(locus: Locus, r1: str, r2: Optional[str], name: str, outdir: str, window: int,
                  platform: str, n_processes: int = 1, min_avg_q: int = 0, exclude_bp: int = 20,
                  extra: Optional[List[str]] = None, count_substitutions: bool = False) -> List[str]:
    nw = NW_PRESETS[platform]
    cmd = ["CRISPResso", "--fastq_r1", r1]
    if r2:
        cmd += ["--fastq_r2", r2]
    cmd += ["--amplicon_seq", locus.amplicon, "--amplicon_name", locus.name,
            "--guide_seq", locus.guide,
            "--quantification_window_size", str(window),
            "--quantification_window_center", str(locus.cleavage_offset),
            "--exclude_bp_from_left", str(exclude_bp), "--exclude_bp_from_right", str(exclude_bp),
            *([] if count_substitutions else ["--ignore_substitutions"]),
            "--needleman_wunsch_gap_open", str(nw["gap_open"]),
            "--needleman_wunsch_gap_extend", str(nw["gap_extend"]),
            "--needleman_wunsch_gap_incentive", str(nw["gap_incentive"]),
            "--name", name, "--output_folder", outdir, "-p", str(n_processes),
            "--place_report_in_output_folder"]
    if locus.coding_seq:
        cmd += ["--coding_seq", locus.coding_seq]
    if min_avg_q:
        cmd += ["--min_average_read_quality", str(min_avg_q)]
    if extra:
        cmd += extra
    return cmd


def run(cmd: List[str], log: Path) -> Path:
    """Run CRISPResso2; on a process-pool permission failure retry with -p 1."""
    def _go(c):
        with open(log, "w") as fh:
            return subprocess.run(c, stdout=fh, stderr=subprocess.STDOUT).returncode
    rc = _go(cmd)
    name = cmd[cmd.index("--name") + 1]
    out = Path(cmd[cmd.index("--output_folder") + 1]) / f"CRISPResso_on_{name}"
    alle = out / "Alleles_frequency_table.zip"
    if (rc != 0 or not alle.exists()) and "SC_SEM_NSEMS_MAX" in log.read_text(errors="ignore"):
        i = cmd.index("-p")
        cmd = cmd[:i + 1] + ["1"] + cmd[i + 2:]
        rc = _go(cmd)
    if not alle.exists():
        tail = "\n".join(log.read_text(errors="ignore").splitlines()[-25:])
        raise RuntimeError(f"CRISPResso2 failed (exit {rc}); see {log}\n{tail}")
    return out


def read_quant(out: Path) -> Dict[str, str]:
    f = out / "CRISPResso_quantification_of_editing_frequency.txt"
    if not f.exists():
        return {}
    rows = list(csv.DictReader(open(f), delimiter="\t"))
    return rows[0] if rows else {}


def read_mapping(out: Path) -> Dict[str, str]:
    f = out / "CRISPResso_mapping_statistics.txt"
    if not f.exists():
        return {}
    rows = list(csv.DictReader(open(f), delimiter="\t"))
    return rows[0] if rows else {}


def iter_alleles(out: Path):
    with zipfile.ZipFile(out / "Alleles_frequency_table.zip") as z:
        name = [n for n in z.namelist() if n.endswith(".txt")][0]
        with z.open(name) as fh:
            for row in csv.DictReader(io.TextIOWrapper(fh), delimiter="\t"):
                yield (row["Aligned_Sequence"], row["Reference_Sequence"], int(row["#Reads"]),
                       int(float(row.get("n_mutated", 0) or 0)))


def rescore(out: Path, locus: Locus, windows=(1, 5, 10, 0), flank: int = 15, exclude_bp: int = 20,
            count_substitutions: bool = False):
    """Re-score the allele table. Returns ({w: Counter}, allele_rows) where
    allele_rows lists the per-allele calls at the first window.  With
    count_substitutions, reads whose only change is a substitution inside the
    CRISPResso2 quantification window (n_mutated > 0; first window only) are
    counted as modified, category 'substitution'."""
    coding = (locus.coding_start, locus.coding_end) if locus.coding_seq else None
    res = {w: Counter() for w in windows}
    allele_rows = []
    L, cut = locus.length, locus.cut
    for aln, ref, n, n_mut in iter_alleles(out):
        indels, rs, re_ = parse_alignment(ref, aln)
        for w in windows:
            c = res[w]
            c["aligned"] += n
            if w > 0:
                ok = rs <= cut - w - flank and re_ >= cut + w + flank
            else:
                ok = rs <= exclude_bp and re_ >= L - exclude_bp
            if not ok:
                c["not_covering_window"] += n
                continue
            call = classify(indels, locus.amplicon, cut, w, coding, exclude_bp, exclude_bp)
            if count_substitutions and w == windows[0] and not call.modified and n_mut > 0:
                call = type(call)(True, "substitution", 0, 0, "sub")
            c["informative"] += n
            c[call.category] += n
            if w == windows[0]:
                allele_rows.append((call.alleles or "WT", call.category, call.net, n, aln))
    return res, allele_rows
