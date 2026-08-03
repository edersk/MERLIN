# =============================================================================
# vmatchPattern2.R - Indel-aware vectorized pattern matching
# =============================================================================

#' Vectorized pattern matching with indel support
#'
#' Finds all occurrences of a pattern in a set of subject sequences,
#' supporting mismatches and indels. Extended version of
#' \code{Biostrings::vmatchPattern} provided by Herve Pages.
#'
#' @param pattern A \code{\link[Biostrings]{DNAString}} pattern to search for.
#' @param subject An \code{\link[Biostrings]{XStringSet}} or character string.
#' @param max.mismatch Maximum number of mismatching letters allowed (default: 0).
#' @param min.mismatch Minimum number of mismatching letters allowed (default: 0).
#' @param with.indels If TRUE, indels are allowed (default: FALSE).
#' @param fixed If TRUE, IUPAC ambiguity codes match only themselves (default: TRUE).
#' @param algorithm Algorithm to use: "auto", "naive-exact", "naive-inexact",
#'   "boyer-moore", "shift-or", or "indels" (default: "auto").
#'
#' @return A list of \code{\link[IRanges]{IRanges}} objects, one per subject sequence,
#'   indicating match positions.
#'
#' Originally provided by Herve Pages on
#' \url{https://support.bioconductor.org/p/58350/}.
#'
#' @examples
#' library(Biostrings)
#' x <- DNAStringSet(c("AAGCGCGATATG", "GCNNNATCCCC"))
#' vmatchPattern2(DNAString("GCNNNAT"), x, max.mismatch = 1)
#'
#' @export
vmatchPattern2 <- function(pattern, subject,
                           max.mismatch = 0, min.mismatch = 0,
                           with.indels = FALSE, fixed = TRUE,
                           algorithm = "auto") {
    # Coerce a plain character subject into the XStringSet the C routine expects.
    if (!methods::is(subject, "XStringSet")) {
        subject <- Biostrings:::XStringSet(NULL, subject)
    }

    # Normalise every argument through the same internal helpers that
    # Biostrings::vmatchPattern() uses, so behaviour matches the public function.
    algo <- Biostrings:::normargAlgorithm(algorithm)
    if (Biostrings:::isCharacterAlgo(algo)) {
        stop(
            "'subject' must be a single (non-empty) string ",
            "for this algorithm"
        )
    }
    pattern <- Biostrings:::normargPattern(pattern, subject)
    max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
    min.mismatch <- Biostrings:::normargMinMismatch(
        min.mismatch,
        max.mismatch
    )
    with.indels <- Biostrings:::normargWithIndels(with.indels)
    fixed <- Biostrings:::normargFixed(fixed, subject)
    algo <- Biostrings:::selectAlgo(
        algo, pattern,
        max.mismatch, min.mismatch,
        with.indels, fixed
    )

    # Run the vectorized match in C; results come back as parallel start/width
    # lists (one entry per subject sequence).
    C_ans <- .Call2("XStringSet_vmatch_pattern", pattern, subject,
        max.mismatch, min.mismatch,
        with.indels, fixed, algo,
        "MATCHES_AS_RANGES",
        PACKAGE = "Biostrings"
    )

    # Flatten to a single IRanges, then re-split it to mirror the subject layout.
    unlisted_ans <- IRanges::IRanges(
        start = unlist(C_ans[[1L]], use.names = FALSE),
        width = unlist(C_ans[[2L]], use.names = FALSE)
    )
    # S4Vectors::relist is not exported; fetch it from the namespace directly to
    # avoid a ::: call while still using the method that understands C_ans.
    relist <- get("relist", envir = asNamespace("S4Vectors"))
    relist(unlisted_ans, C_ans[[1L]])
}
