# =============================================================================
# zzz.R - Namespace imports and package load hook for MERLIN
# =============================================================================
#
# The @importFrom directives below make every external function used by the
# pipeline resolvable from within the package namespace, so the function bodies
# can call them unqualified. A handful of calls are qualified inline instead:
# Biostrings ::: internals in vmatchPattern2(), pwalign::pairwiseAlignment() /
# pwalign::score(), parallel::detectCores(), and S4Vectors' relist().
# =============================================================================

#' @importFrom Biostrings DNAString DNAStringSet alphabetFrequency extractAt
#'   nucleotideSubstitutionMatrix quality reverseComplement subseq
#' @importFrom BiocGenerics start end width
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors .Call2
#' @importFrom ShortRead FastqStreamer ShortReadQ id sread yield
#' @importFrom data.table data.table fread fwrite copy setcolorder setkey
#'   setnames setorder uniqueN fcase fifelse is.data.table := .GRP .I .N .SD
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics abline axis grid hist legend lines par plot plot.new
#'   polygon segments text
#' @importFrom inflection uik
#' @importFrom methods is
#' @importFrom parallel detectCores mclapply
#' @importFrom pwalign pairwiseAlignment score
#' @importFrom stats density loess median predict quantile setNames smooth.spline
#' @importFrom stringdist stringdistmatrix
#' @importFrom stringr str_split
NULL

# Verify that the non-exported Biostrings helpers relied on by vmatchPattern2()
# are present in the installed Biostrings. They are reached with ::: because no
# public API exposes indel-aware vectorized matching; warn early if the internal
# names have changed rather than failing deep inside the pipeline.
.onLoad <- function(libname, pkgname) {
    required_internals <- c(
        "XStringSet", "normargAlgorithm", "isCharacterAlgo",
        "normargPattern", "normargMaxMismatch", "normargMinMismatch",
        "normargWithIndels", "normargFixed", "selectAlgo"
    )
    for (fn in required_internals) {
        if (!exists(fn, envir = asNamespace("Biostrings"))) {
            warning(sprintf(
                "Biostrings internal '%s' not found. vmatchPattern2() may not work.",
                fn
            ), call. = FALSE)
        }
    }
}
