# crispr-indel

Measure CRISPR knockout efficiency from amplicon sequencing: the percentage of
reads carrying an indel at the Cas9 cut site, and the percentage that are
frameshifted. Built around [CRISPResso2](https://github.com/pinellolab/CRISPResso2),
with **Oxford Nanopore single-end reads as the primary input** and Illumina
paired-end supported.

Preconfigured for the human **PROCR** exon-3 guide
(`config/PROCR_g1.yaml`); any other locus needs one YAML file or three
command-line arguments.

```
                 % modified   % frameshift   control-corrected
 unedited WT          0.79           0.79                0.00
 PROCR KO pool       58.28          48.10               57.95
```

---

## Quick start

```bash
git clone https://github.com/pmtran5884/crispr-indel.git
cd crispr-indel

# 1. install (one command, ~2 min; works on Linux and Apple Silicon)
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate crispr-indel

# 2. check it works on the bundled synthetic data
./run_example.sh

# 3. run your data
python run_crispr_indel.py \
    --locus config/PROCR_g1.yaml \
    --control WT.fastq.gz \
    KO_clone1.fastq.gz KO_clone2.fastq.gz \
    -o results/
```

Open `results/report.html`. The key numbers are in `results/editing_summary.tsv`.

A whole directory of per-barcode FASTQs works too:

```bash
python run_crispr_indel.py --locus config/PROCR_g1.yaml \
    --control fastq_pass/barcode01 fastq_pass/barcode02 fastq_pass/barcode03 -o results/
```

No conda? `pip install .` then `crispr-indel ...`; the tool runs with its
built-in aligner and says so if CRISPResso2 is absent.
Docker: `docker build -t crispr-indel . && docker run --rm -v "$PWD":/data -w /data crispr-indel --help`.

---

## What you need

| | |
|---|---|
| **Reads** | one FASTQ per sample (ONT) or R1/R2 (Illumina), gzipped or not |
| **A control** | an unedited sample (parental line / non-targeting guide). Strongly recommended for ONT — see below |
| **A locus file** | `config/PROCR_g1.yaml` is included; for another target see *Another guide or gene* |

Nothing else — no reference genome, no BAM, no manual trimming. Barcodes,
adapters and read orientation are handled internally.

### Read-length requirement

Reads must physically span the cut site. The PROCR amplicon is 428 bp with the
cut 232 bp from the forward primer, so:

* **ONT:** fine — reads cover the full amplicon.
* **Illumina:** 2×250 (mates merge), **2×150 will not work** and the tool stops
  with an explanation rather than reporting a wrong number.

---

## Why an unedited control matters (nanopore)

Nanopore basecalling makes indel errors in homopolymers, and the PROCR cut site
sits between two `TTT` runs (`CCTTTG|TGAGTTT`). An unedited sample therefore
shows a non-zero "editing" background:

| chemistry / basecalling | apparent editing in an unedited control |
|---|---|
| R10.4.1, SUP (recommended) | ~0.8–1% |
| R9.4.1 or fast basecalling | ~3% |

The control sample removes this exactly: `(treated − control) / (100 − control)`.
Without a control, subtract nothing and treat anything below ~3% as noise.

---

## Output

```
results/
├── editing_summary.tsv          <- the headline table, one row per sample
├── all_metrics.tsv              <- both engines x windows (+/-1, 5, 10, whole amplicon)
├── report.html                  <- self-contained: figures, top alleles, notes
├── outcomes.png                 <- stacked frameshift / in-frame bar chart
├── <sample>_indel_sizes.png     <- indel size distribution
├── run.log
└── <sample>/
    ├── top_alleles.tsv          <- most frequent alleles with categories
    ├── amplicon_reads.fastq.gz  <- oriented, primer-trimmed reads (ONT)
    ├── builtin_read_calls.tsv   <- per-read call from the built-in engine
    └── CRISPResso_on_<sample>/  <- the full CRISPResso2 output and its own report
```

Columns in `editing_summary.tsv`:

| column | meaning |
|---|---|
| `informative_reads` | reads that actually span the cut-site window (the denominator) |
| `pct_modified` | % of those with an indel at the cut site |
| `pct_modified_ci95_lo/hi` | Wilson 95% CI (read sampling only) |
| `pct_frameshift` | % with a coding indel not divisible by 3 — the loss-of-function proxy |
| `pct_inframe` | % with an in-frame coding indel (may retain function) |
| `pct_modified_corrected` | **the number to report** when a control was given |
| `builtin_pct_modified` | the independent engine's value; a >2-point difference triggers a warning |

---

## Interpreting the numbers

* **"Editing efficiency"** in most papers = `pct_modified` (control-corrected).
* **Knockout efficiency** is better approximated by `pct_frameshift`; in-frame
  indels may leave a functional protein, so report both.
* These are **allele** frequencies, not cell frequencies. A polyclonal pool at
  60% modified alleles is not 60% knockout cells.
* High-efficiency samples read **slightly low** on nanopore (about 3% relative
  on R10, 9% on R9) because a sequencing error at the cut can mask a real
  indel. Treat >80% values as lower bounds.
* Alleles with large deletions that remove a primer site are not amplified and
  are invisible to this (or any) amplicon assay.

---

## Validation

`tools/validate.py` simulates reads with a known editing rate (realistic SpCas9
outcome spectrum, nanopore homopolymer/truncation/chimera errors), runs the
full pipeline, and compares with per-read truth.

![validation](docs/validation.png)

| true % edited | measured, R10.4.1-like | measured, R9.4.1-like |
|---|---|---|
| 0 | 0.00 | 0.00 |
| ~10 | 9.50 | 8.37 |
| ~25 | 24.18 | 23.72 |
| ~50 | 48.73 | 45.01 |
| ~75 | 72.07 | 68.35 |
| ~90 | 86.39 | 81.99 |

Reproduce with `python tools/validate.py --out validation_run` (~3 min).
`pytest -q` runs 11 unit/integration tests including these recovery checks.

---

## Another guide or gene

Either write a YAML like `config/PROCR_g1.yaml`, or let the tool build one by
in-silico PCR against Ensembl GRCh38:

```bash
python tools/build_amplicon.py --gene PROCR \
    --forward TTCCCTGACCTAGACTCCCC --reverse CCCCGTCCCATCATACCTTT \
    --guide CTCCGGCCGGAAACTCACAA --name PROCR_g1 --out config/PROCR_g1.yaml
```

It finds the amplicon, places the cut site, checks the PAM, and extracts the
coding sequence of the exon containing the cut for frameshift calls.

For a one-off run without a file:

```bash
python run_crispr_indel.py --amplicon <AMPLICON_SEQ> --guide <20nt_PROTOSPACER> \
    --control WT.fastq.gz KO.fastq.gz -o results/
```

---

## Useful options

| option | default | when to change |
|---|---|---|
| `-w, --window` | 1 | `-w 5` or `-w 10` to count indels further from the cut |
| `--min-q` | 10 | raise to 12–15 to be stricter with ONT reads |
| `--count-substitutions` | off | base editors, or to recover ONT reads where an error masks an indel |
| `--engine` | both | `quick` to skip CRISPResso2, `crispresso` to skip the cross-check |
| `--max-reads` | all | e.g. `--max-reads 20000` for a fast first look |
| `-p, --threads` | 1 | CRISPResso2 processes |
| `--platform` | auto | force `ont` / `illumina` if detection by read length is wrong |

---

## The preconfigured PROCR locus

| | |
|---|---|
| amplicon | chr20:35,176,035–35,176,462 (GRCh38), 428 bp |
| primers | F `TTCCCTGACCTAGACTCCCC` / R `CCCCGTCCCATCATACCTTT` |
| guide | `CTCCGGCCGGAAACTCACAA`, minus strand, PAM `AGG` |
| cut site | chr20:35,176,266^35,176,267 — 232 bp from F, 196 bp from R |
| context | `...GAGCTCCTTTG | TGAGTTTCCGGCCGG...` |
| coding | exon 3 of ENST00000216968 (PROCR canonical), cut is in-frame-sensitive coding sequence |

Verified against Ensembl REST: the primer pair and guide each occur exactly
once in the region, and the guide sits 232 bp into the amplicon — close to
centred, as intended.

---

## Method

See [docs/METHOD.md](docs/METHOD.md) for the algorithm, parameter choices
(including why gap penalties differ for nanopore), the control-correction
model, validation details, and limitations.

## Citing

If you use this, cite CRISPResso2 (Clement et al., *Nat Biotechnol* 2019) and,
for nanopore amplicon analysis, McFarlane et al., *BMC Res Notes* 2024.

## License

MIT
