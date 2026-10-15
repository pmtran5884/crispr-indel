"""Independent indel quantification engine (cross-check / fallback).

Aligns every read to the amplicon with an affine-gap semi-global alignment
(parasail; the read is aligned end-to-end, amplicon ends are free), extracts
indels, and classifies the read with the same repeat-aware window rules as the
CRISPResso2 allele-table parser (crispr_indel.indels).

Single reads (ONT after crispr_indel.prep, or Illumina single-end): a read is
informative if it covers the window +/- `flank` bp.
Paired Illumina reads (not merged): each mate is aligned separately; if both
mates are informative and disagree, the mate with the higher mean base quality
over the window wins (and the pair is counted as discordant).
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Iterable, Iterator, List, Optional, Tuple

import parasail

from .indels import ReadCall, classify, parse_alignment
from .locus import Locus, revcomp
from .preflight import primer_orientation, read_fastq

# Scoring presets.  Illumina: CRISPResso2-like (EDNAFULL 5/-4, gap open 20, extend 2).
# ONT: minimap2 map-ont-like (2/-4, gap 4+2k).  With Illumina-style penalties a
# nanopore sequencing indel a few bp from a real Cas9 indel is often "cancelled"
# into a run of mismatches (two gaps cost more than 3-4 mismatches), hiding the
# edit; cheap gaps relative to mismatches keep both indels.
SCORING = {
    "illumina": dict(match=5, mismatch=-4, gap_open=20, gap_extend=2, min_score_per_base=2.0),
    "ont": dict(match=2, mismatch=-4, gap_open=6, gap_extend=2, min_score_per_base=0.6),
}
_MATS = {}


def _matrix(match, mismatch):
    k = (match, mismatch)
    if k not in _MATS:
        _MATS[k] = parasail.matrix_create("ACGTN", match, mismatch)
    return _MATS[k]
ADAPTERS = ("AGATCGGAAGAGC", "CTGTCTCTTATACACATCT")


@dataclass
class Aln:
    ok: bool
    indels: list = None
    ref_start: int = 0
    ref_end: int = 0
    score: int = 0
    qual_window: float = 0.0


def align(seq: str, qual: str, locus: Locus, oriented: bool = False, platform: str = "illumina") -> Aln:
    sc = SCORING[platform]
    M = _matrix(sc["match"], sc["mismatch"])
    ref = locus.amplicon
    for ad in ADAPTERS:
        i = seq.find(ad)
        if i >= 0:
            seq, qual = seq[:i], qual[:i]
    if len(seq) < 30:
        return Aln(False)
    cands = [(seq, qual)]
    if not oriented:
        o = primer_orientation(seq, locus)
        if o == "R":
            cands = [(revcomp(seq), qual[::-1])]
        elif o == "":
            cands = [(seq, qual), (revcomp(seq), qual[::-1])]
    best = None
    for s, q in cands:
        r = parasail.sg_dx_trace_scan_16(s, ref, sc["gap_open"], sc["gap_extend"], M)
        if best is None or r.score > best[0].score:
            best = (r, s, q)
    r, s, q = best
    if r.score < sc["min_score_per_base"] * len(s):
        return Aln(False, score=r.score)
    tb = r.traceback
    indels, rs, re_ = parse_alignment(tb.ref, tb.query)
    # mean read quality over the bases aligned to the cut +/- 10 bp
    qw, rp, qp = [], 0, 0
    for a, b in zip(tb.ref, tb.query):
        if a != "-" and b != "-" and abs(rp - locus.cut) <= 10 and qp < len(q):
            qw.append(ord(q[qp]) - 33)
        rp += a != "-"
        qp += b != "-"
    return Aln(True, indels, rs, re_, r.score, sum(qw) / len(qw) if qw else 0.0)


def informative(a: Aln, locus: Locus, w: int, flank: int) -> bool:
    if not a.ok:
        return False
    if w <= 0:
        return a.ref_start <= 20 and a.ref_end >= locus.length - 20
    return a.ref_start <= locus.cut - w - flank and a.ref_end >= locus.cut + w + flank


@dataclass
class QuickResult:
    counts: Counter
    calls: List[Tuple[str, str, ReadCall]]  # (read name, status, call)


def _call(a: Aln, locus: Locus, w: int, excl: Tuple[int, int]) -> ReadCall:
    coding = (locus.coding_start, locus.coding_end) if locus.coding_seq else None
    return classify(a.indels, locus.amplicon, locus.cut, w, coding, excl[0], excl[1])


def run_single(fastq: str, locus: Locus, w: int = 1, flank: int = 15, oriented: bool = False,
               excl: Tuple[int, int] = (20, 20), max_reads: Optional[int] = None,
               keep_calls: bool = True, platform: str = "ont") -> QuickResult:
    counts: Counter = Counter()
    calls = []
    for i, (name, s, q) in enumerate(read_fastq(fastq)):
        if max_reads and i >= max_reads:
            break
        counts["reads"] += 1
        a = align(s, q, locus, oriented=oriented, platform=platform)
        if not a.ok:
            counts["unaligned"] += 1
            if keep_calls:
                calls.append((name, "unaligned", None))
            continue
        if not informative(a, locus, w, flank):
            counts["not_covering_window"] += 1
            if keep_calls:
                calls.append((name, "not_covering", None))
            continue
        c = _call(a, locus, w, excl)
        counts["informative"] += 1
        counts[c.category] += 1
        if keep_calls:
            calls.append((name, "ok", c))
    return QuickResult(counts, calls)


def run_paired(r1: str, r2: str, locus: Locus, w: int = 1, flank: int = 15,
               excl: Tuple[int, int] = (20, 20), max_reads: Optional[int] = None,
               keep_calls: bool = True) -> QuickResult:
    counts: Counter = Counter()
    calls = []
    for i, ((n1, s1, q1), (n2, s2, q2)) in enumerate(zip(read_fastq(r1), read_fastq(r2))):
        if max_reads and i >= max_reads:
            break
        counts["reads"] += 1
        a1, a2 = align(s1, q1, locus), align(s2, q2, locus)
        i1, i2 = informative(a1, locus, w, flank), informative(a2, locus, w, flank)
        if not (i1 or i2):
            counts["unaligned" if not (a1.ok or a2.ok) else "not_covering_window"] += 1
            if keep_calls:
                calls.append((n1, "uninformative", None))
            continue
        c1 = _call(a1, locus, w, excl) if i1 else None
        c2 = _call(a2, locus, w, excl) if i2 else None
        if c1 and c2:
            if c1.modified != c2.modified:
                counts["discordant_pairs"] += 1
            c = c1 if a1.qual_window >= a2.qual_window else c2
        else:
            c = c1 or c2
        counts["informative"] += 1
        counts[c.category] += 1
        if keep_calls:
            calls.append((n1, "ok", c))
    return QuickResult(counts, calls)


def summarize_counts(c: Counter) -> dict:
    n = c.get("informative", 0)
    mod = n - c.get("unmodified", 0)
    pct = (lambda x: 100.0 * x / n) if n else (lambda x: float("nan"))
    return {
        "reads": c.get("reads", 0),
        "informative_reads": n,
        "pct_modified": pct(mod),
        "pct_frameshift": pct(c.get("frameshift", 0)),
        "pct_inframe": pct(c.get("in-frame", 0)),
        "pct_exon_boundary": pct(c.get("exon-boundary", 0)),
        "pct_noncoding": pct(c.get("noncoding", 0)),
        "unaligned": c.get("unaligned", 0),
        "not_covering_window": c.get("not_covering_window", 0),
        "discordant_pairs": c.get("discordant_pairs", 0),
    }
