"""Long-read (ONT/PacBio) preprocessing: primer-anchored extraction.

Nanopore amplicon reads carry barcode/adapter flanks of variable length, come
in both orientations, and include truncated molecules and head-to-tail
chimeras.  CRISPResso2 assumes reads are adapter-trimmed and will otherwise
score truncated reads as 'unmodified', biasing editing downwards.

For every read we locate the forward primer and the reverse-complemented
reverse primer by semi-global alignment (tolerant of ONT errors), orient the
read to the amplicon (+) strand, and keep only the primer-to-primer product.
Reads are binned as:
  ok           both primers found, product length plausible, mean Q >= min_q
  one_primer   truncated molecule / failed read (cannot see whole amplicon)
  no_primer    off-target product, primer dimer, junk
  too_long     product > amplicon + max_insertion (concatemer / chimera)
  low_quality  mean Q of product < min_q
"""
from __future__ import annotations

import gzip
from collections import Counter
from dataclasses import dataclass
from typing import Optional, Tuple

import parasail

from .locus import Locus, revcomp
from .preflight import read_fastq

_PM = parasail.matrix_create("ACGTN", 2, -3)


def find_primer(seq: str, primer: str, lo: int, hi: int) -> Tuple[int, int, int]:
    """Best placement of `primer` (fully aligned) inside seq[lo:hi].
    Returns (score, start, end) in seq coordinates (end exclusive)."""
    lo, hi = max(0, lo), min(len(seq), hi)
    if hi - lo < len(primer) // 2:
        return -999, -1, -1
    r = parasail.sg_dx_trace_scan_16(primer, seq[lo:hi], 5, 2, _PM)
    q, ref = r.traceback.query, r.traceback.ref
    # columns where the primer has a base
    cols = [i for i, c in enumerate(q) if c != "-"]
    if not cols:
        return -999, -1, -1
    first, last = cols[0], cols[-1]
    start = lo + sum(1 for c in ref[:first] if c != "-")
    end = lo + sum(1 for c in ref[:last + 1] if c != "-")
    return r.score, start, end


@dataclass
class PrepResult:
    status: str
    seq: str = ""
    qual: str = ""
    strand: str = ""


def extract(seq: str, qual: str, locus: Locus, search: int = 250, min_primer_frac: float = 0.7,
            min_q: float = 0.0, max_insertion: int = 150) -> PrepResult:
    F, R = locus.forward_primer, locus.reverse_primer
    rcF, rcR = revcomp(F), revcomp(R)
    thrF, thrR = min_primer_frac * 2 * len(F), min_primer_frac * 2 * len(R)
    n = len(seq)
    best = None
    for strand, p5, p3, t5, t3 in (("+", F, rcR, thrF, thrR), ("-", R, rcF, thrR, thrF)):
        s5, a5, _ = find_primer(seq, p5, 0, search + len(p5))
        s3, _, b3 = find_primer(seq, p3, n - search - len(p3), n)
        hit5, hit3 = s5 >= t5, s3 >= t3
        if hit5 and hit3 and b3 - a5 >= len(F) + len(R):
            cand = (2, s5 + s3, strand, a5, b3)
        elif hit5 or hit3:
            cand = (1, max(s5, s3), strand, -1, -1)
        else:
            cand = (0, max(s5, s3), strand, -1, -1)
        if best is None or cand[:2] > best[:2]:
            best = cand
    level, _, strand, a, b = best
    if level == 0:
        return PrepResult("no_primer")
    if level == 1:
        return PrepResult("one_primer")
    prod, pq = seq[a:b], qual[a:b]
    if strand == "-":
        prod, pq = revcomp(prod), pq[::-1]
    if len(prod) > locus.length + max_insertion:
        return PrepResult("too_long")
    if min_q and pq:
        mq = sum(ord(c) - 33 for c in pq) / len(pq)
        if mq < min_q:
            return PrepResult("low_quality")
    return PrepResult("ok", prod, pq, strand)


def prep_fastq(inp: str, out: str, locus: Locus, min_q: float = 0.0, search: int = 250,
               min_primer_frac: float = 0.7, max_insertion: int = 150,
               max_reads: Optional[int] = None) -> Counter:
    """Write oriented, primer-trimmed reads to `out` (gz). Returns status counts
    (plus strand counts under keys 'strand+' / 'strand-')."""
    counts: Counter = Counter()
    with gzip.open(out, "wt") as fo:
        for i, (name, s, q) in enumerate(read_fastq(inp)):
            if max_reads and i >= max_reads:
                break
            r = extract(s, q, locus, search, min_primer_frac, min_q, max_insertion)
            counts[r.status] += 1
            if r.status == "ok":
                counts["strand" + r.strand] += 1
                fo.write(f"@{name}\n{r.seq}\n+\n{r.qual}\n")
    return counts
