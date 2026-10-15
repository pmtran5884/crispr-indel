"""FASTQ sniffing and pre-flight checks.

The most common way an amplicon-seq indel analysis goes wrong is not the
algorithm but the reads: they are too short to reach the cut site, the mates do
not overlap enough to merge, or the library is non-directional.  These checks
run in seconds on the first few thousand reads and either pass, warn, or stop
with an explanation.
"""
from __future__ import annotations

import gzip
import statistics
from dataclasses import dataclass, field
from itertools import islice
from typing import Iterator, Optional, Tuple

from .locus import Locus, revcomp


def open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def read_fastq(path: str) -> Iterator[Tuple[str, str, str]]:
    """Yield (name, seq, qual). Minimal, fast, no dependencies."""
    with open_text(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().rstrip()
            fh.readline()
            q = fh.readline().rstrip()
            if not q:
                raise ValueError(f"truncated FASTQ record in {path}")
            yield h[1:].split()[0], s, q


def _hamming_prefix(read: str, primer: str) -> int:
    n = len(primer)
    return sum(a != b for a, b in zip(read[:n], primer)) + max(0, n - len(read))


def primer_orientation(seq: str, locus: Locus, max_mm: int = 3) -> str:
    """'F' if the read starts with the forward primer, 'R' if with the reverse
    primer, '' otherwise.  Allows a few mismatches and a 0-3 nt offset (stagger/
    leftover adapter)."""
    best = ("", max_mm + 1)
    for off in range(0, 4):
        s = seq[off:]
        for tag, p in (("F", locus.forward_primer), ("R", locus.reverse_primer)):
            if p:
                d = _hamming_prefix(s, p)
                if d < best[1]:
                    best = (tag, d)
    return best[0]


@dataclass
class FastqStats:
    path: str
    n_sampled: int = 0
    median_len: float = 0
    p90_len: float = 0
    max_len: int = 0
    mean_q: float = 0
    frac_F: float = 0
    frac_R: float = 0

    def line(self) -> str:
        return (f"{self.path}: n={self.n_sampled} median_len={self.median_len:.0f} max_len={self.max_len} "
                f"mean_Q={self.mean_q:.1f} starts_with_F={100*self.frac_F:.1f}% starts_with_R={100*self.frac_R:.1f}%")


def sniff(path: str, locus: Locus, n: int = 5000) -> FastqStats:
    lens, qs, f, r = [], [], 0, 0
    for _, s, q in islice(read_fastq(path), n):
        lens.append(len(s))
        if q:
            qs.append(sum(ord(c) - 33 for c in q) / len(q))
        o = primer_orientation(s, locus)
        f += o == "F"
        r += o == "R"
    st = FastqStats(path=path, n_sampled=len(lens))
    if lens:
        sl = sorted(lens)
        st.median_len = statistics.median(sl)
        st.p90_len = sl[int(0.9 * (len(sl) - 1))]
        st.max_len = sl[-1]
        st.mean_q = statistics.fmean(qs) if qs else 0
        st.frac_F, st.frac_R = f / len(lens), r / len(lens)
    return st


@dataclass
class Preflight:
    mode: str                 # 'merged' | 'single' | 'fail'
    messages: list = field(default_factory=list)
    stats: list = field(default_factory=list)
    ok: bool = True


def preflight(locus: Locus, r1: str, r2: Optional[str], window: int, flank: int = 15,
              min_overlap: int = 10, n: int = 5000) -> Preflight:
    """Decide how reads can be used to quantify editing at the cut site.

    A read starting at the forward primer must extend at least cut + window +
    flank bases to see the whole quantification window with anchoring sequence
    on the far side; a read starting at the reverse primer needs
    (L - cut) + window + flank.  Paired reads can be merged if R1 + R2 >= L +
    min_overlap.
    """
    L, cut = locus.length, locus.cut
    need_from_F = cut + window + flank
    need_from_R = (L - cut) + window + flank
    s1 = sniff(r1, locus, n)
    stats = [s1]
    s2 = sniff(r2, locus, n) if r2 else None
    if s2:
        stats.append(s2)
    msgs = [f"Locus: {locus.summary()}"] + [s.line() for s in stats]
    msgs.append(f"A read must extend >= {need_from_F} nt from the F primer or >= {need_from_R} nt from the "
                f"R primer to cover the cut site +/- {window} bp with {flank} bp anchor.")
    oriented = s1.frac_F + s1.frac_R
    if s1.n_sampled == 0:
        return Preflight("fail", msgs + ["ERROR: R1 is empty."], stats, False)
    if oriented < 0.5:
        msgs.append(f"WARNING: only {100*oriented:.0f}% of R1 reads start with a primer. Reads may still contain "
                    "adapters/barcodes, or this is not the expected amplicon (check primers / sample).")
    nondirectional = s1.frac_F > 0.1 and s1.frac_R > 0.1
    if nondirectional:
        msgs.append(f"Library is non-directional (R1 starts with F primer in {100*s1.frac_F:.0f}% and with R primer "
                    f"in {100*s1.frac_R:.0f}% of reads) - typical of ligation-based amplicon services; handled.")
    l1 = s1.median_len
    if s2:
        l2 = s2.median_len
        if l1 + l2 >= L + min_overlap:
            msgs.append(f"Paired-end {l1:.0f}+{l2:.0f} nt vs {L} bp amplicon: expected mate overlap "
                        f"{l1 + l2 - L:.0f} nt -> reads will be MERGED (fastp) before alignment.")
            if l1 + l2 - L < 30:
                msgs.append("WARNING: overlap < 30 nt; insertions longer than ~"
                            f"{int(l1 + l2 - L - min_overlap)} nt will prevent merging and be under-counted.")
            return Preflight("merged", msgs, stats, True)
        msgs.append(f"Mates ({l1:.0f}+{l2:.0f} nt) are too short to overlap across a {L} bp amplicon.")
    # single-read usage
    covers_F = l1 >= need_from_F
    covers_R = l1 >= need_from_R
    if (nondirectional and covers_F and covers_R) or (not nondirectional and (
            (s1.frac_F >= s1.frac_R and covers_F) or (s1.frac_R > s1.frac_F and covers_R))):
        msgs.append("Single reads reach across the cut site -> running in SINGLE-READ mode on R1.")
        return Preflight("single", msgs, stats, True)
    msgs.append(f"ERROR: reads cannot cover the cut site. Cut is {cut} bp from the F primer and {L - cut} bp "
                f"from the R primer, but reads are ~{l1:.0f} nt"
                + (f" and mates cannot be merged" if s2 else "") +
                ". Use longer reads (e.g. 2x250 for this amplicon) or a shorter amplicon.")
    if nondirectional and (covers_F or covers_R):
        msgs.append("Only reads from one primer end reach the cut site; you may force this with --mode single "
                    "(reads from the other end will be unaligned/uninformative and lower the read count).")
    return Preflight("fail", msgs, stats, False)
