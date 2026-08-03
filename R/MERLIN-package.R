#' MERLIN: Mutation-Enriched RNA Profiling via Long-Read Integration
#'
#' MERLIN genotypes single nucleotide variants (SNVs) from Oxford Nanopore
#' long-read single-cell RNA sequencing data. It streams FASTQ files, extracts
#' cell barcodes and UMIs, calls the base at one or two user-defined regions of
#' interest (ROIs), corrects UMIs, and classifies each cell as mutant or
#' wild-type with optional SNP-based phasing.
#'
#' @section Pipeline:
#' The package exposes three functions that are run in sequence:
#' \enumerate{
#'   \item \code{\link{detect_snv}} -- stream FASTQ, detect the Read 1 adapter,
#'     extract the cell barcode (CBC) and UMI, require a poly(T) tail, and call
#'     the base at the major (and optional minor) ROI.
#'   \item \code{\link{summarize_snv}} -- correct UMIs by Levenshtein distance,
#'     aggregate reads per cell barcode, detect knee/inflection thresholds, and
#'     report SNP allele frequencies to assess phasing feasibility.
#'   \item \code{\link{flag_snv}} -- classify each cell as \code{MUT}, \code{WT},
#'     or \code{MUT_WT}, with an optional SNP-phased call that separates
#'     definitive from ambiguous wild-type cells.
#' }
#' The indel-aware matcher \code{\link{vmatchPattern2}} underpins all sequence
#' searches in the pipeline.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

# Register the non-standard-evaluation (data.table) column names used across the
# pipeline so that R CMD check does not flag them as undefined global variables.
#' @importFrom utils globalVariables
globalVariables(c(
    ".", ".GRP", ".I", ".N", ".SD", "..cols_keep",
    # detect_snv() working columns
    "index", "id", "width", "read_id",
    "read1_fwd_end", "read1_rev_end", "read1_end_combined",
    "is_valid_read1", "is_dual_read", "strand", "polyN_exists",
    "CBC", "UMI", "major_base", "minor_base", "major_status", "minor_status",
    # summarize_snv() working columns
    "initial_UMI", "corrected_UMI", "corrected", "N", "n_reads", "rank",
    "n_mut", "n_wt", "n_total", "n_minor_mut", "n_minor_wt",
    "n_major", "n_minor", "major_class", "minor_threshold", "minor_class",
    "y_val", "group", "n_umis", "n_cells_with_status", "status", "roi",
    "pct_reads", "pct_umis",
    # summary / flag_snv() output columns
    "cbc", "umis", "n_umi", "n_umi_mut", "n_umi_wt",
    "n_snp", "n_snp_mut", "n_snp_wt",
    "has_np_data", "mutation_call", "phased_call",
    "major_call_simple", "phased_call_simple", "barcode_10x"
))

# Signal that this package is aware of data.table's reference semantics.
.datatable.aware <- TRUE
