"""Indel extraction and classification from a pairwise alignment.

Shared by the CRISPResso2 allele-table parser and the independent quick
aligner, so both engines are scored by exactly the same rules:

* An indel is *in the window* if ANY equivalent placement of it (shifting it
  left/right through repeated sequence, e.g. a 1-nt deletion in a TTT run)
  overlaps the quantification window around the cut site.  This mirrors what
  CRISPResso2's gap incentive does at alignment time and makes the call
  independent of how an aligner happens to place a gap in a repeat.
* Window half-size w: bases amplicon[cut-w : cut+w] (CRISPResso2 convention;
  w=1 -> the two bases flanking the cut). Insertions count if they can be
  placed strictly inside the window (w=1 -> exactly at the cut).
* w = 0 means "anywhere in the amplicon except the excluded primer ends".
* Substitutions are ignored for indel calls (they are dominated by PCR /
  sequencing error in Cas9 nuclease experiments).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple


@dataclass(frozen=True)
class Indel:
    kind: str      # 'D' or 'I'
    pos: int       # D: 0-based ref start; I: insertion point (between pos-1 and pos)
    length: int
    seq: str = ""  # inserted bases (I) or deleted ref bases (D)

    def label(self) -> str:
        return f"{self.pos}{self.kind}{self.length}" + (f"{self.seq}" if self.kind == "I" else "")


def parse_alignment(ref_aln: str, read_aln: str) -> Tuple[List[Indel], int, int]:
    """Return (indels, ref_start, ref_end) from gapped alignment strings.

    Terminal gaps are NOT indels: leading/trailing '-' in the read mean the read
    does not cover that part of the reference; leading/trailing '-' in the
    reference are read overhang.  ref_start/ref_end (0-based, half-open) is the
    reference span covered by the read.
    """
    assert len(ref_aln) == len(read_aln)
    n = len(ref_aln)
    # trim terminal columns where either side is a gap at the ends
    i, j = 0, n
    while i < j and (read_aln[i] == "-" or ref_aln[i] == "-"):
        i += 1
    while j > i and (read_aln[j - 1] == "-" or ref_aln[j - 1] == "-"):
        j -= 1
    ref_pos = sum(1 for c in ref_aln[:i] if c != "-")
    ref_start = ref_pos
    indels: List[Indel] = []
    k = i
    while k < j:
        a, b = ref_aln[k], read_aln[k]
        if a != "-" and b != "-":
            ref_pos += 1
            k += 1
        elif b == "-":  # deletion from read
            s = ref_pos
            seq = []
            while k < j and read_aln[k] == "-" and ref_aln[k] != "-":
                seq.append(ref_aln[k]); ref_pos += 1; k += 1
            indels.append(Indel("D", s, len(seq), "".join(seq)))
        else:  # insertion in read
            seq = []
            while k < j and ref_aln[k] == "-":
                seq.append(read_aln[k]); k += 1
            indels.append(Indel("I", ref_pos, len(seq), "".join(seq)))
    return indels, ref_start, ref_pos


def equivalent_range(ind: Indel, ref: str) -> Tuple[int, int]:
    """Leftmost and rightmost equivalent positions of an indel in `ref`."""
    L = len(ref)
    if ind.kind == "D":
        s, e = ind.pos, ind.pos + ind.length
        sl, el = s, e
        while sl > 0 and ref[sl - 1] == ref[el - 1]:
            sl -= 1; el -= 1
        sr, er = s, e
        while er < L and ref[sr] == ref[er]:
            sr += 1; er += 1
        return sl, sr
    # insertion: rotate inserted sequence through flanking repeat
    x, p = ind.seq, ind.pos
    xl, pl = x, p
    while pl > 0 and xl and ref[pl - 1] == xl[-1]:
        xl = xl[-1] + xl[:-1]; pl -= 1
    xr, pr = x, p
    while pr < L and xr and ref[pr] == xr[0]:
        xr = xr[1:] + xr[0]; pr += 1
    return pl, pr


def in_window(ind: Indel, ref: str, cut: int, w: int, excl_left: int = 0, excl_right: int = 0) -> bool:
    lo, hi = equivalent_range(ind, ref)
    if w <= 0:
        a, b = excl_left, len(ref) - excl_right
        if ind.kind == "D":
            return lo < b and hi + ind.length > a
        return hi >= a and lo <= b
    if ind.kind == "D":
        return lo < cut + w and hi + ind.length > cut - w
    return lo <= cut + w - 1 and hi >= cut - w + 1


@dataclass
class ReadCall:
    modified: bool
    category: str          # 'unmodified' | 'frameshift' | 'in-frame' | 'exon-boundary' | 'noncoding'
    net: int               # net length change of in-window indels (ins - del)
    max_del: int
    alleles: str           # compact allele label, e.g. '232I1T;225D7'


def classify(indels: List[Indel], ref: str, cut: int, w: int,
             coding: Optional[Tuple[int, int]] = None,
             excl_left: int = 0, excl_right: int = 0, large_del: int = 50) -> ReadCall:
    hits = [d for d in indels if in_window(d, ref, cut, w, excl_left, excl_right)]
    if not hits:
        return ReadCall(False, "unmodified", 0, 0, "")
    net = sum(d.length if d.kind == "I" else -d.length for d in hits)
    max_del = max([d.length for d in hits if d.kind == "D"], default=0)
    label = ";".join(d.label() for d in hits)
    if coding is None or coding[0] < 0:
        return ReadCall(True, "noncoding", net, max_del, label)
    cs, ce = coding
    coding_hits, boundary = [], False
    for d in hits:
        if d.kind == "D":
            s, e = d.pos, d.pos + d.length
            if s < ce and e > cs:
                if s < cs or e > ce:
                    boundary = True
                coding_hits.append(d)
        elif cs < d.pos < ce:
            coding_hits.append(d)
    if boundary:
        cat = "exon-boundary"
    elif not coding_hits:
        cat = "noncoding"
    else:
        cnet = sum(d.length if d.kind == "I" else -d.length for d in coding_hits)
        cat = "frameshift" if cnet % 3 else "in-frame"
    return ReadCall(True, cat, net, max_del, label)
