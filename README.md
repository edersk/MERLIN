# MERLIN

**M**utation-**E**nriched **R**NA profiling via **L**ong-read **In**tegration

MERLIN genotypes single nucleotide variants (SNVs) from Oxford Nanopore
long-read single-cell RNA sequencing data. It streams FASTQ files to extract
cell barcodes and UMIs, calls the base at user-defined regions of interest,
performs UMI error correction, and classifies cells as mutant or wild-type with
optional SNP-based phasing.

## Installation

```r
# Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "IRanges", "S4Vectors", "ShortRead",
                       "BiocGenerics", "pwalign"))

# Install MERLIN from GitHub
# install.packages("remotes")
remotes::install_github("edersk/MERLIN")
```

MERLIN also depends on the CRAN packages `data.table`, `stringr`, `stringdist`
and `inflection`, which are installed automatically.

## Pipeline

Three functions are run in sequence:

```
FASTQ --> detect_snv() --> summarize_snv() --> flag_snv()
```

1. **`detect_snv()`** — streams the FASTQ, detects the Read 1 adapter, extracts
   the cell barcode (CBC) and UMI, requires a poly(T) sequnce, and calls the base at
   the major region of interest (ROI) and an optional minor ROI (a SNP for phasing).
2. **`summarize_snv()`** — corrects UMIs by Levenshtein distance, aggregates
   reads per cell, detects knee/inflection thresholds, and reports SNP allele
   frequencies to gauge phasing feasibility.
3. **`flag_snv()`** — classifies each cell as `MUT`, `WT`, or `MUT_WT`, with an
   optional SNP-phased call that separates definitive from ambiguous wild-type.

## Quick start

Point the pipeline at your own FASTQ and cell-barcode whitelist:

```r
library(MERLIN)

fastq    <- "path/to/reads.fastq.gz"   # Nanopore long reads
barcodes <- "path/to/barcodes.csv"     # cell-barcode whitelist (10x -1 suffix OK)
outdir   <- tempfile("merlin_"); dir.create(outdir)

# Step 1: detect SNVs from FASTQ
detect_snv(
    session_name       = "experiment_name",
    path_input_file    = "path/to/experiment.fastq.gz",
    path_output_folder = "outdir/",
    read1_seq          = "CTACACGACGCTCTTCCGATCT", # adapter sequence
    read1_mismatch     = 3L,
    major_roi_seq      = "GAGATTTCNCTGTAGCT",   # 'N' marks the SNV position
    major_roi_mut      = "T",
    major_roi_wt       = "A",
    major_roi_mismatch = 3L,
    minor_roi_seq      = "AGAACAATNCCAAATGC",   # optional SNP for phasing
    minor_roi_mut      = "C",
    minor_roi_wt       = "T",
    minor_roi_mismatch = 3L,
    n_cores            = 4L,
    n_reads_per_chunk = 1e5    
)

# Step 2: correct UMIs and summarize per cell
res <- summarize_snv(
    session_name       = "experiment_name",
    path_snv_table     = "outdir/experiment_name_snv_table.csv"),
    path_barcodes      = "path/to/barcodes.csv"     # 10x or individual experiment cell-barcode whitelist
    path_output_folder = "output/"
    umi_mismatch = 2L,
    n_cores = 4L,
)

res$knee_points        # suggested read-per-UMI thresholds
res$phasing_message    # whether phasing looks viable for this dataset

# Step 3: classify cells (optionally phased)
flags <- flag_snv(
    path_snv_summary = "output/experiment_name_snv_summary.csv",
    path_barcode_filter = NULL, # optinal CBC filter `path/to/barcodes.csv`
    path_output_folder = "output/",
    session_name = "experiment_name",
    threshold = 10,    # possible recommended thresholds shown summarize_snv
    use_phasing = FALSE,
    minor_on_major = "MUT"    # Phasing SNP status on major SNV allele
)

table(flags$mutation_call)
table(flags$phased_call)
```

## Input requirements

- **FASTQ** (`.fastq` or `.fastq.gz`) of Oxford Nanopore long reads.
- **ROI sequence** — the region flanking the variant, with a single `N` at the
  SNV position, e.g. `GAGATTTCNCTGTAGCT`. 
- **Barcode whitelist** — a `.csv` (or plain-text) list of expected cell
  barcodes; only the first 16 bp are used for matching.
- A **minor ROI** is optional; provide one only if you want SNP-based phasing.

## Key parameters

More customizable variables like cell barcode and UMI length; polyN base, length, mismatch rate; consensus threshold; etc. are provided in the detailed documentation of the functions.

| Parameter                       | Function        | Meaning                                                                                      |
| ------------------------------- | --------------- | -------------------------------------------------------------------------------------------- |
| `read1_seq`, `read1_mismatch`   | `detect_snv`    | Adapter sequence preceding the barcode and its mismatch tolerance (number of bases)          |
| `read1_within_range`            | `detect_snv`    | How far into the read to search for the adapter (default: 70)                                |
| `major_roi_seq`                 | `detect_snv`    | ROI containing a single `N` at the SNV position                                              |
| `major_roi_mut`, `major_roi_wt` | `detect_snv`    | Bases labelling mutant vs. wild-type                                                         |
| `minor_roi_seq`                 | `detect_snv`    | Optional SNP sequence used for phasing containing a single `N` at the SNV position           |
| `polyN_base`, `polyN_length`    | `detect_snv`    | Tail required after the UMI (default poly-T, length at least 10)                             |
| `len_CBC`, `len_UMI`            | `detect_snv`    | Cell barcode and UMI lengths (default: 16 / 12)                                              |
| `n_reads_per_chunk`, `n_cores`  | `detect_snv`    | Streaming chunk size and number of cores                                                     |
| `umi_mismatch`                  | `summarize_snv` | Max edit distance for UMI collapsing                                                         |
| `threshold`                     | `flag_snv`      | Minimum total reads per UMI to keep                                                          |
| `consensus_threshold`           | `flag_snv`      | Fraction a UMI must lean to be called (default: 0.8)                                         |
| `use_phasing`, `minor_on_major` | `flag_snv`      | Use the heterozygous minor SNP; which minor SNP status (MUT or WT) sits on the mutant allele |

## Outputs

| File                             | Written by      | Contents                                                                                                                                |
| -------------------------------- | --------------- | --------------------------------------------------------------------------------------------------------------------------------------- |
| `<session>_snv_table.csv`        | `detect_snv`    | One row per read: `read_id`, `CBC`, `UMI`, `strand`, `major_base`/`major_status`, `minor_base`/`minor_status`                           |
| `<session>_snv_summary.csv`      | `summarize_snv` | One row per cell, `;`-separated: `cbc;umis;n_umi;n_umi_mut;n_umi_wt;n_snp;n_snp_mut;n_snp_wt` (count columns are `,`-separated per UMI) |
| `<session>_umi_frequencies.csv`  | `summarize_snv` | Ranked reads-per-UMI, for choosing a threshold                                                                                          |
| `<session>_*_roi_frequency.csv`  | `summarize_snv` | Major/minor allele frequency tables                                                                                                     |
| `<session>_umi_distribution.pdf` | `summarize_snv` | Diagnostic plots (rank-abundance, histogram, cumulative, violins)                                                                       |
| `<session>_flags.csv`            | `flag_snv`      | Per-cell genotype: `cbc`, `barcode_10x`, `has_np_data`, `mutation_call`, `major_call_simple`, `phased_call`, `phased_call_simple`       |

## Genotype categories

**`mutation_call`** (major ROI): `MUT`, `WT`, `MUT_WT` (heterozygous — both
alleles), `below_threshold` (no UMI cleared the read threshold), or `no_np_data`
(barcode present in the whitelist but had no Nanopore reads). `major_call_simple`
collapses these to `MUT` / `WT`.

**`phased_call`** (only meaningful with `use_phasing = TRUE`): the minor SNP is
used to qualify the major call

| `phased_call`          | `phased_call_simple` | Meaning                                                  |
| ---------------------- | -------------------- | -------------------------------------------------------- |
| `MUT_mnsnp`            | `MUT`                | Mutant, confirmed by the expected minor SNP              |
| `MUT_non_mnsnp`        | `MUT`                | Mutant, minor SNP not observed                           |
| `MUT_both`             | `MUT`                | Mutant with a mix of observed and non observed minor SNP |
| `def_WT`               | `defWT`              | Definitive wild-type (phasing-consistent)                |
| `amb_WT`               | `ambWT`              | Ambiguous wild-type (could not be phased)                |
| `MUT_ERROR` / `ERROR1` | —                    | Minor SNP contradicts the major call                     |

## Choosing a threshold

`summarize_snv()` returns `res$knee_points` and writes
`<session>_umi_distribution.pdf`. Pick the `flag_snv()` `threshold` within the
suggested range: a **higher** threshold is more specific (fewer, higher-
confidence calls); a **lower** threshold is more sensitivea with possible false-positive calls. 

## Notes

- `vmatchPattern2()` relies on non-exported Biostrings C helpers (via `:::`),
  which produces NOTEs in `R CMD check`. 

## License

MIT © 2026 Sebastian Eder
