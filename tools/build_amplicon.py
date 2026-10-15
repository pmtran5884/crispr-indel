#!/usr/bin/env python
"""Build a locus YAML from primers + guide + gene symbol (in-silico PCR on GRCh38).

Fetches the gene region (+/- 5 kb) from the Ensembl REST API, locates the
primer pair, extracts the amplicon in forward-primer -> reverse-primer
orientation, validates the guide/PAM/cut site, and takes the coding part of the
canonical-transcript exon that contains the cut (for frameshift calls).

  python tools/build_amplicon.py --gene PROCR --forward TTCCCTGACCTAGACTCCCC \\
      --reverse CCCCGTCCCATCATACCTTT --guide CTCCGGCCGGAAACTCACAA --name PROCR_g1 \\
      --out config/PROCR_g1.yaml

Needs internet access to rest.ensembl.org; species defaults to human.
"""
import argparse
import datetime
import json
import os
import re
import sys
import urllib.request

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from crispr_indel.locus import Locus, revcomp  # noqa: E402

API = "https://rest.ensembl.org"


def get(path, ctype="application/json"):
    req = urllib.request.Request(API + path, headers={"Content-Type": ctype})
    with urllib.request.urlopen(req, timeout=60) as r:
        data = r.read().decode()
    return json.loads(data) if ctype == "application/json" else data.strip()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--species", default="homo_sapiens")
    ap.add_argument("--forward", required=True)
    ap.add_argument("--reverse", required=True)
    ap.add_argument("--guide", required=True)
    ap.add_argument("--name", default=None)
    ap.add_argument("--pad", type=int, default=5000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    F, R, G = (x.strip().upper() for x in (a.forward, a.reverse, a.guide))
    g = get(f"/lookup/symbol/{a.species}/{a.gene}?expand=1")
    chrom, s0, e0 = g["seq_region_name"], g["start"] - a.pad, g["end"] + a.pad
    seq = get(f"/sequence/region/{a.species}/{chrom}:{s0}..{e0}:1", "text/plain").upper()
    hits = []
    for strand, f, r in (("+", F, revcomp(R)), ("-", R, revcomp(F))):
        for m in re.finditer(f"(?={f})", seq):
            j = seq.find(r, m.start() + len(f))
            if 0 <= j < m.start() + 5000:
                hits.append((strand, m.start(), j + len(r)))
    if len(hits) != 1:
        sys.exit(f"expected exactly one amplicon, found {len(hits)}: {hits[:5]}")
    strand, i, j = hits[0]
    amp = seq[i:j] if strand == "+" else revcomp(seq[i:j])
    gstart, gend = s0 + i, s0 + j - 1
    loc = Locus(name=a.name or f"{a.gene}_amplicon", amplicon=amp, guide=G, forward_primer=F,
                reverse_primer=R).validate()
    # genomic coordinate of each amplicon index
    def gpos(k):
        return gstart + k if strand == "+" else gend - k
    cut_g = sorted([gpos(loc.cut - 1), gpos(loc.cut)])
    tx = [t for t in g.get("Transcript", []) if t.get("is_canonical")]
    coding, tinfo = "", {}
    if tx and "Translation" in tx[0]:
        t = tx[0]
        cds_s, cds_e = t["Translation"]["start"], t["Translation"]["end"]
        for n, ex in enumerate(sorted(t["Exon"], key=lambda e: e["start"]), 1):
            cs, ce = max(ex["start"], cds_s, gstart), min(ex["end"], cds_e, gend)
            if cs <= cut_g[0] <= ce:
                if strand == "+":
                    coding = amp[cs - gstart: ce - gstart + 1]
                else:
                    coding = amp[gend - ce: gend - cs + 1]
                rank = n if t["strand"] == 1 else len(t["Exon"]) - n + 1
                tinfo = dict(transcript=t["id"], exon=rank)
    if coding:
        loc.coding_seq = coding
        loc.validate()
    lines = [f"# {a.gene} locus built by tools/build_amplicon.py from Ensembl REST ({datetime.date.today()})",
             f"name: {loc.name}", f"gene: {a.gene}", f"genome: {g.get('assembly_name', '')}",
             f"amplicon_coords: \"chr{chrom}:{gstart}-{gend}\"   # amplicon given in F->R orientation "
             f"({'+' if strand == '+' else '-'} genomic strand)",
             f"forward_primer: {F}", f"reverse_primer: {R}", f"guide: {G}",
             f"pam: NGG", f"pam_observed: {loc.meta['pam_observed']}",
             f"cut_site_genomic: \"chr{chrom}:{cut_g[0]}^{cut_g[1]}\"",
             f"cut_index_in_amplicon: {loc.cut}", f"amplicon_length: {loc.length}",
             f"amplicon: {amp}"]
    if tinfo:
        lines += [f"transcript: {tinfo['transcript']}", f"exon: {tinfo['exon']}"]
    lines += [f"coding_seq: {coding}" if coding else "# coding_seq: (cut not in a coding exon)"]
    open(a.out, "w").write("\n".join(lines) + "\n")
    print(loc.summary())
    for w in loc.warnings:
        print("WARNING:", w)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
