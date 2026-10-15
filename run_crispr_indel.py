#!/usr/bin/env python
"""Run without installing:  python run_crispr_indel.py --locus config/PROCR_g1.yaml --control WT.fastq.gz KO.fastq.gz"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from crispr_indel.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
