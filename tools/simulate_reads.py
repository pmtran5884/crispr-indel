#!/usr/bin/env python
"""Simulate paired-end amplicon reads with a known CRISPR editing rate.

Used for validation/tests (no real data is shipped).  The edited allele pool
follows a typical SpCas9 NHEJ/MMEJ outcome spectrum: +1 templated insertions
dominate, followed by -1, short deletions, microhomology deletions, a few
larger insertions, large (50-90 bp) deletions and complex indels.  Reads are
non-directional (R1 starts at either primer, as with ligation-based amplicon
services), carry a position-dependent substitution error profile, rare
sequencing indels, and ~1% junk reads.

Writes <prefix>_R1.fastq.gz, <prefix>_R2.fastq.gz and <prefix>_truth.tsv
(one row per read pair: allele label, whether it is edited, true category).

Example:
  python tools/simulate_reads.py --locus config/PROCR_g1.yaml --edit 0.5 \
      --pairs 20000 --read-len 250 --out example_data/sim50
"""
import argparse
import gzip
import math
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from crispr_indel.indels import Indel, classify  # noqa: E402
from crispr_indel.locus import load_locus, revcomp  # noqa: E402

ADAPT_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPT_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"


def apply(ref, indels):
    out, p = [], 0
    for d in sorted(indels, key=lambda x: (x.pos, x.kind == "D")):  # insertion before a deletion at the same pos
        out.append(ref[p:d.pos])
        if d.kind == "D":
            p = d.pos + d.length
        else:
            out.append(d.seq); p = d.pos
    out.append(ref[p:])
    return "".join(out)


def allele_pool(ref, cut, rng, coding):
    """Return list of (indels, weight) for edited alleles."""
    pool = []
    B = "ACGT"
    # +1 templated insertion (duplicate base 5' of the cut on the guide strand) + other +1
    pool.append(([Indel("I", cut, 1, ref[cut - 1])], 0.25))
    pool.append(([Indel("I", cut, 1, ref[cut])], 0.08))
    for b in B:
        if b not in (ref[cut - 1], ref[cut]):
            pool.append(([Indel("I", cut, 1, b)], 0.02))
    # -1
    pool.append(([Indel("D", cut - 1, 1, ref[cut - 1])], 0.08))
    pool.append(([Indel("D", cut, 1, ref[cut])], 0.06))
    # short deletions 2-20 spanning the cut
    for _ in range(14):
        n = min(20, 2 + int(rng.expovariate(1 / 5)))
        s = rng.randint(cut - n + 1, cut)  # deletion [s, s+n) contains cut-1 or cut
        s = max(s, cut - n + 1)
        pool.append(([Indel("D", s, n, ref[s:s + n])], 0.20 / 14))
    # microhomology deletions: repeats of length >=2 on both sides of the cut
    mh = []
    for k in range(2, 6):
        for i in range(cut - 25, cut):
            for j in range(cut, cut + 25):
                if ref[i:i + k] == ref[j:j + k] and j - i >= 3 and i < cut <= j:
                    mh.append(Indel("D", i, j - i, ref[i:j]))
    rng.shuffle(mh)
    for d in mh[:6]:
        pool.append(([d], 0.10 / max(1, min(6, len(mh)))))
    # 2-10 bp insertions
    for _ in range(4):
        n = rng.randint(2, 10)
        pool.append(([Indel("I", cut, n, "".join(rng.choice(B) for _ in range(n)))], 0.04 / 4))
    # large deletions 50-90 bp that keep both primers and stay in the exon
    for _ in range(4):
        n = rng.randint(50, 90)
        lo = max(coding[0] + 1, cut - n + 1) if coding else cut - n + 1
        hi = min(cut, (coding[1] - n - 1) if coding else cut)
        s = rng.randint(lo, max(lo, hi))
        pool.append(([Indel("D", s, n, ref[s:s + n])], 0.06 / 4))
    # complex: deletion + insertion at the cut
    for _ in range(2):
        n = rng.randint(3, 12); s = cut - rng.randint(1, n - 1)
        ins = "".join(rng.choice(B) for _ in range(rng.randint(1, 4)))
        pool.append(([Indel("D", s, n, ref[s:s + n]), Indel("I", s, len(ins), ins)], 0.03 / 2))
    tot = sum(w for _, w in pool)
    return [(a, w / tot) for a, w in pool]


def err_rate(i, L, scale):
    return scale * (0.001 + 0.012 * (i / L) ** 2)


def sequence(frag, read_len, adapter, rng, scale, indel_rate=2e-5):
    s = (frag + adapter + "G" * read_len)[:read_len]
    out, qs = [], []
    for i, c in enumerate(s):
        if rng.random() < indel_rate:
            continue  # sequencing deletion
        e = err_rate(i, read_len, scale)
        if rng.random() < e:
            c = rng.choice([b for b in "ACGT" if b != c])
        q = int(max(2, min(41, -10 * math.log10(e) + rng.gauss(0, 2))))
        out.append(c); qs.append(chr(q + 33))
    return "".join(out), "".join(qs)


def homopolymer_runs(seq):
    """Return list of (start, length) runs."""
    runs, i = [], 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        runs.append((i, j - i)); i = j
    return runs


def sequence_ont(frag, rng, scale):
    """ONT-like read: random barcode/adapter flanks, substitution + indel errors,
    with indel errors concentrated in homopolymers (length-dependent), and
    per-base qualities ~Q10-30.  scale=1 approximates R10.4.1 + SUP basecalling
    (~1% total error); scale=4 approximates R9.4.1 / fast basecalling."""
    left = "".join(rng.choice("ACGT") for _ in range(rng.randint(40, 90)))
    right = "".join(rng.choice("ACGT") for _ in range(rng.randint(30, 70)))
    out, qs = [], []
    p_sub, p_ins, p_del = 0.004 * scale, 0.002 * scale, 0.003 * scale
    for start, n in homopolymer_runs(left + frag + right):
        b = (left + frag + right)[start]
        # homopolymer length error: P grows with run length
        p_hp = scale * {1: 0.0, 2: 0.004, 3: 0.012, 4: 0.03}.get(n, 0.06)
        r = rng.random()
        if r < p_hp * 0.65:
            n2 = n - 1
        elif r < p_hp:
            n2 = n + 1
        else:
            n2 = n
        for _ in range(n2):
            r = rng.random()
            if r < p_del:
                continue
            c = b
            if r < p_del + p_sub:
                c = rng.choice([x for x in "ACGT" if x != b])
            out.append(c)
            q = int(max(3, min(50, rng.gauss(24, 7) if c == b else rng.gauss(10, 3))))
            qs.append(chr(q + 33))
            if rng.random() < p_ins:
                out.append(rng.choice("ACGT")); qs.append(chr(int(max(3, rng.gauss(9, 3))) + 33))
    return "".join(out), "".join(qs)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--edit", type=float, required=True, help="fraction of read pairs from edited alleles")
    ap.add_argument("--pairs", type=int, default=20000)
    ap.add_argument("--read-len", type=int, default=250)
    ap.add_argument("--error-scale", type=float, default=1.0)
    ap.add_argument("--junk", type=float, default=0.01)
    ap.add_argument("--directional", action="store_true", help="R1 always starts at the forward primer")
    ap.add_argument("--single-end", action="store_true")
    ap.add_argument("--platform", choices=["illumina", "ont"], default="illumina",
                    help="ont: single-end full-length reads with an ONT-like error model")
    ap.add_argument("--truncated", type=float, default=0.03, help="ONT: fraction of reads truncated mid-amplicon")
    ap.add_argument("--chimera", type=float, default=0.01, help="ONT: fraction of head-to-tail chimeric reads")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True, help="output prefix")
    a = ap.parse_args(argv)
    rng = random.Random(a.seed)
    loc = load_locus(a.locus)
    ref, cut = loc.amplicon, loc.cut
    coding = (loc.coding_start, loc.coding_end) if loc.coding_seq else None
    pool = allele_pool(ref, cut, random.Random(a.seed + 1000), coding)
    alleles, weights = zip(*pool)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    f1 = gzip.open(a.out + "_R1.fastq.gz", "wt")
    ont = a.platform == "ont"
    f2 = None if (a.single_end or ont) else gzip.open(a.out + "_R2.fastq.gz", "wt")
    tr = open(a.out + "_truth.tsv", "w")
    tr.write("read\tedited\tcategory\tallele\n")
    for n in range(a.pairs):
        name = f"sim{n}"
        if rng.random() < a.junk:
            frag = "".join(rng.choice("ACGT") for _ in range(rng.randint(60, 300)))
            lab, edited, cat = "junk", "NA", "junk"
        elif rng.random() < a.edit:
            ind = rng.choices(alleles, weights)[0]
            frag = apply(ref, ind)
            rc_ = classify(list(ind), ref, cut, 1, coding)
            lab, edited, cat = ";".join(d.label() for d in ind), "1", rc_.category
        else:
            frag, lab, edited, cat = ref, "WT", "0", "unmodified"
        if not a.directional and rng.random() < 0.5:
            frag = revcomp(frag)
        if ont:
            if cat != "junk" and rng.random() < a.truncated:
                frag = frag[:rng.randint(80, len(frag) - 60)]
                cat, edited = "truncated", "NA"
            elif cat != "junk" and rng.random() < a.chimera:
                frag = frag + revcomp(ref)
                cat, edited = "chimera", "NA"
            s1, q1 = sequence_ont(frag, rng, a.error_scale)
        else:
            s1, q1 = sequence(frag, a.read_len, ADAPT_R1, rng, a.error_scale)
        f1.write(f"@{name} 1:N:0:1\n{s1}\n+\n{q1}\n")
        if f2:
            s2, q2 = sequence(revcomp(frag), a.read_len, ADAPT_R2, rng, a.error_scale * 1.5)
            f2.write(f"@{name} 2:N:0:1\n{s2}\n+\n{q2}\n")
        tr.write(f"{name}\t{edited}\t{cat}\t{lab}\n")
    f1.close(); tr.close()
    if f2:
        f2.close()
    print(f"wrote {a.pairs} reads to {a.out}_R1.fastq.gz" + ("" if f2 is None else " (+R2)"))


if __name__ == "__main__":
    main()
