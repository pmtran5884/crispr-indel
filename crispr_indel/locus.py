"""Locus definition: amplicon, primers, guide, cut site, coding sequence.

A locus can be given as a YAML file (see config/PROCR_g1.yaml) or built on the
fly from an amplicon sequence + guide.  Everything is validated against the
amplicon so that a typo in a guide or primer fails loudly instead of silently
producing a wrong cut site.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]


def _clean(s: Optional[str]) -> str:
    return re.sub(r"\s+", "", s or "").upper()


@dataclass
class Locus:
    name: str
    amplicon: str
    guide: str
    forward_primer: str = ""
    reverse_primer: str = ""
    coding_seq: str = ""
    pam: str = "NGG"
    cleavage_offset: int = -3  # SpCas9: cut 3 nt 5' of the PAM (CRISPResso2 convention)
    meta: dict = field(default_factory=dict)

    # filled in by validate()
    guide_strand: str = ""
    guide_start: int = -1   # 0-based start of the protospacer in amplicon (+ strand coords)
    cut: int = -1           # 0-based: cut lies between amplicon[cut-1] and amplicon[cut]
    coding_start: int = -1  # 0-based half-open interval of coding_seq in amplicon
    coding_end: int = -1
    warnings: list = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.amplicon)

    def validate(self) -> "Locus":
        self.amplicon = _clean(self.amplicon)
        self.guide = _clean(self.guide)
        self.forward_primer = _clean(self.forward_primer)
        self.reverse_primer = _clean(self.reverse_primer)
        self.coding_seq = _clean(self.coding_seq)
        if not re.fullmatch(r"[ACGTN]+", self.amplicon):
            raise ValueError("amplicon contains non-ACGTN characters")
        g = self.guide
        fwd = [m.start() for m in re.finditer(f"(?={g})", self.amplicon)]
        rev = [m.start() for m in re.finditer(f"(?={revcomp(g)})", self.amplicon)]
        if len(fwd) + len(rev) == 0:
            raise ValueError("guide not found in amplicon on either strand")
        if len(fwd) + len(rev) > 1:
            raise ValueError("guide matches the amplicon more than once; cut site ambiguous")
        pam_len = len(self.pam)
        pam_re = self.pam.replace("N", "[ACGT]").replace("R", "[AG]")
        if fwd:
            p = fwd[0]
            self.guide_strand, self.guide_start = "+", p
            pam_seq = self.amplicon[p + len(g): p + len(g) + pam_len]
            self.cut = p + len(g) + self.cleavage_offset
        else:
            p = rev[0]
            self.guide_strand, self.guide_start = "-", p
            pam_seq = revcomp(self.amplicon[max(0, p - pam_len): p])
            self.cut = p - self.cleavage_offset
        self.meta["pam_observed"] = pam_seq
        if not re.fullmatch(pam_re, pam_seq):
            self.warnings.append(f"sequence adjacent to guide is '{pam_seq}', not a {self.pam} PAM")
        if self.forward_primer and not self.amplicon.startswith(self.forward_primer):
            if self.forward_primer in self.amplicon:
                self.warnings.append("forward primer is internal to the amplicon (amplicon should start at the primer)")
            else:
                raise ValueError("forward primer not found at the start of the amplicon")
        if self.reverse_primer and not self.amplicon.endswith(revcomp(self.reverse_primer)):
            if revcomp(self.reverse_primer) in self.amplicon:
                self.warnings.append("reverse primer is internal to the amplicon (amplicon should end at its reverse complement)")
            else:
                raise ValueError("reverse complement of the reverse primer not found at the end of the amplicon")
        if self.coding_seq:
            i = self.amplicon.find(self.coding_seq)
            if i < 0:
                raise ValueError("coding_seq is not a substring of the amplicon (+ strand)")
            self.coding_start, self.coding_end = i, i + len(self.coding_seq)
            if not (self.coding_start <= self.cut <= self.coding_end):
                self.warnings.append("cut site is outside coding_seq; frameshift calls will mostly be 'noncoding'")
        return self

    def summary(self) -> str:
        a = self.amplicon
        return (f"{self.name}: amplicon {self.length} bp; guide on {self.guide_strand} strand, "
                f"PAM {self.meta.get('pam_observed')}; cut between amplicon[{self.cut - 1}] and "
                f"amplicon[{self.cut}] ({a[self.cut - 5:self.cut]}|{a[self.cut:self.cut + 5]}); "
                f"{self.cut} bp from the forward end, {self.length - self.cut} bp from the reverse end"
                + (f"; coding {self.coding_start}-{self.coding_end}" if self.coding_seq else ""))


def load_locus(path: Optional[str] = None, amplicon: Optional[str] = None,
               guide: Optional[str] = None, coding_seq: Optional[str] = None,
               forward_primer: Optional[str] = None, reverse_primer: Optional[str] = None,
               name: Optional[str] = None) -> Locus:
    """Load a locus from YAML, with optional command-line overrides."""
    d: dict = {}
    if path:
        import yaml
        d = yaml.safe_load(Path(path).read_text()) or {}
    known = {"name", "amplicon", "guide", "forward_primer", "reverse_primer", "coding_seq",
             "pam", "cleavage_offset"}
    kw = {k: d[k] for k in known if k in d and d[k] is not None}
    for k, v in dict(amplicon=amplicon, guide=guide, coding_seq=coding_seq,
                     forward_primer=forward_primer, reverse_primer=reverse_primer, name=name).items():
        if v:
            kw[k] = v
    if "amplicon" not in kw or "guide" not in kw:
        raise ValueError("need at least an amplicon and a guide (via --locus YAML or --amplicon/--guide)")
    kw.setdefault("name", "amplicon")
    meta = {k: v for k, v in d.items() if k not in known}
    return Locus(meta=meta, **kw).validate()
