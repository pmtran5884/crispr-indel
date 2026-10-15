"""Tests: locus validation, indel normalisation, ONT preprocessing, recovery of
known editing rates by both engines, and the read-length guard.

Run:  pytest -q          (CRISPResso2 tests are skipped if it is not installed)
"""
import csv
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from crispr_indel import crispresso as CR  # noqa: E402
from crispr_indel.indels import Indel, classify, equivalent_range, in_window, parse_alignment  # noqa: E402
from crispr_indel.locus import load_locus, revcomp  # noqa: E402
from crispr_indel.prep import extract, prep_fastq  # noqa: E402
from crispr_indel.quick import run_single  # noqa: E402

LOCUS = os.path.join(ROOT, "config", "PROCR_g1.yaml")
SIM = os.path.join(ROOT, "tools", "simulate_reads.py")
CLI = os.path.join(ROOT, "run_crispr_indel.py")


@pytest.fixture(scope="module")
def locus():
    return load_locus(LOCUS)


def sim(tmp, name, edit, n=1500, seed=1, extra=()):
    p = os.path.join(tmp, name)
    subprocess.run([sys.executable, SIM, "--locus", LOCUS, "--edit", str(edit), "--pairs", str(n),
                    "--seed", str(seed), "--out", p, *extra], check=True, stdout=subprocess.DEVNULL)
    rows = [r for r in csv.DictReader(open(p + "_truth.tsv"), delimiter="\t") if r["edited"] in "01"]
    return p, 100 * sum(r["edited"] == "1" for r in rows) / len(rows)


def test_locus_geometry(locus):
    assert locus.length == 428
    assert locus.guide_strand == "-"
    assert locus.meta["pam_observed"] == "AGG"
    assert locus.cut == 232
    assert locus.amplicon.startswith(locus.forward_primer)
    assert locus.amplicon.endswith(revcomp(locus.reverse_primer))
    assert locus.coding_start < locus.cut < locus.coding_end


def test_bad_guide_rejected():
    with pytest.raises(ValueError):
        load_locus(LOCUS, guide="ACGTACGTACGTACGTACGT")


def test_parse_alignment_terminal_gaps_not_indels():
    ref = "ACGTACGTAC"
    ind, s, e = parse_alignment("ACGTACGTAC", "---TAC-TAC")
    assert s == 3 and e == 10
    assert ind == [Indel("D", 6, 1, "G")]


def test_repeat_normalisation(locus):
    ref, cut = locus.amplicon, locus.cut
    # 1-nt deletion inside the TTT run just 5' of the cut: equivalent placements 228..230
    d = Indel("D", 228, 1, "T")
    assert equivalent_range(d, ref) == (228, 230)
    assert not in_window(d, ref, cut, 1)        # does not touch bases 231/232
    assert in_window(d, ref, cut, 5)
    # deleting the G at 231, placed anywhere equivalent, is in the +/-1 window
    assert in_window(Indel("D", 231, 1, "G"), ref, cut, 1)
    # insertion of T at the cut is in the window
    assert in_window(Indel("I", cut, 1, "T"), ref, cut, 1)


def test_frameshift_classification(locus):
    cod = (locus.coding_start, locus.coding_end)
    assert classify([Indel("I", locus.cut, 1, "T")], locus.amplicon, locus.cut, 1, cod).category == "frameshift"
    assert classify([Indel("D", locus.cut - 1, 3, "GTG")], locus.amplicon, locus.cut, 1, cod).category == "in-frame"
    assert classify([], locus.amplicon, locus.cut, 1, cod).category == "unmodified"


def test_ont_prep_bins_reads(tmp_path, locus):
    p, _ = sim(str(tmp_path), "ont", 0.5, extra=("--platform", "ont"))
    truth = {r["read"]: r["category"] for r in csv.DictReader(open(p + "_truth.tsv"), delimiter="\t")}
    from crispr_indel.preflight import read_fastq
    bad_kept = valid_lost = n_valid = 0
    for name, s, q in read_fastq(p + "_R1.fastq.gz"):
        st = extract(s, q, locus).status
        if truth[name] in ("junk", "truncated", "chimera"):
            bad_kept += st == "ok"
        else:
            n_valid += 1
            valid_lost += st != "ok"
    assert bad_kept <= 2
    assert valid_lost / n_valid < 0.02


@pytest.mark.parametrize("edit", [0.0, 0.5, 0.9])
def test_builtin_engine_recovers_ont(tmp_path, locus, edit):
    p, t = sim(str(tmp_path), f"ont{edit}", edit, n=2000, seed=11, extra=("--platform", "ont"))
    prep_fastq(p + "_R1.fastq.gz", str(tmp_path / "prep.fq.gz"), locus)
    c = run_single(str(tmp_path / "prep.fq.gz"), locus, 1, oriented=True, platform="ont").counts
    pct = 100 * (c["informative"] - c["unmodified"]) / c["informative"]
    # R10-like reads: <=2% absolute background, <=~5% relative under-calling at high editing
    assert t - 0.06 * t - 1.0 <= pct <= t + 2.0, (t, pct)


@pytest.mark.skipif(not CR.crispresso_available(), reason="CRISPResso2 not installed")
def test_cli_end_to_end_ont(tmp_path):
    pc, _ = sim(str(tmp_path), "ctrl", 0.0, n=1500, seed=21, extra=("--platform", "ont"))
    pe, t = sim(str(tmp_path), "ko", 0.6, n=1500, seed=22, extra=("--platform", "ont"))
    out = tmp_path / "res"
    subprocess.run([sys.executable, CLI, "--locus", LOCUS, "--control", pc + "_R1.fastq.gz",
                    pe + "_R1.fastq.gz", "-o", str(out)], check=True, stdout=subprocess.DEVNULL)
    rows = {r["sample"]: r for r in csv.DictReader(open(out / "editing_summary.tsv"), delimiter="\t")}
    ko = rows["ko_R1"]
    assert ko["engine"] == "CRISPResso2"
    assert abs(float(ko["pct_modified_corrected"]) - t) < 4.0
    assert abs(float(ko["pct_modified"]) - float(ko["builtin_pct_modified"])) < 2.0
    assert float(rows["ctrl_R1"]["pct_modified"]) < 2.0
    assert (out / "report.html").exists()


def test_short_illumina_reads_fail_preflight(tmp_path):
    p, _ = sim(str(tmp_path), "il150", 0.5, n=500, extra=("--read-len", "150"))
    tsv = tmp_path / "s.tsv"
    tsv.write_text(f"sample\tfastq_r1\tfastq_r2\nx\t{p}_R1.fastq.gz\t{p}_R2.fastq.gz\n")
    r = subprocess.run([sys.executable, CLI, "--locus", LOCUS, "--samples", str(tsv), "-o",
                        str(tmp_path / "o"), "--engine", "quick"], capture_output=True, text=True)
    assert r.returncode != 0
    assert "cannot cover the cut site" in r.stdout + r.stderr
