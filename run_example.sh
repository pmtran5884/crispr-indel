#!/usr/bin/env bash
# Smoke test on the bundled synthetic nanopore data (60.4% true editing, 50.2% frameshift).
set -euo pipefail
cd "$(dirname "$0")"
python run_crispr_indel.py \
  --locus config/PROCR_g1.yaml \
  --control example_data/example_ont_unedited_control.fastq.gz \
  example_data/example_ont_PROCR_KO.fastq.gz \
  -o example_results "$@"
echo; echo "Open example_results/report.html"
