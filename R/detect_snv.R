# =============================================================================
# detect_snv.R - SNV detection from Nanopore single-cell reads (MERLIN)
# =============================================================================
#
# Streams a FASTQ file in chunks and, for each read, locates the Read 1 adapter,
# extracts the cell barcode (CBC) and UMI, requires a poly(T) tail, and calls the
# base found at one or two regions of interest (ROIs). The heavy lifting is done
# with indel-aware vectorized matching (vmatchPattern2) plus a two-pass strategy
# that only falls back to expensive pairwise alignment for reads that do not
# match exactly.
# =============================================================================


# --- Internal helpers --------------------------------------------------------

#' Locate a pattern in every read and return the end position of a unique match
#'
#' Uses indel-aware vectorized matching so all reads are searched in a single
#' call instead of looping. Reads with zero or multiple matches return 0.
#'
#' @param pattern DNAString pattern to search for.
#' @param subject_set DNAStringSet of reads to search within.
#' @param max_mismatch Maximum number of mismatches allowed.
#' @return Integer vector (one per read) with the end coordinate of the unique
#'   match, or 0 where there was not exactly one match.
#' @noRd
.sequence_detection_vectorized <- function(pattern, subject_set, max_mismatch) {
    # Batch match: much faster than per-read mclapply.
    matches <- vmatchPattern2(pattern, subject_set,
        max.mismatch = max_mismatch,
        min.mismatch = 0,
        fixed = TRUE,
        with.indels = TRUE
    )

    # Pre-allocate; only reads with exactly one hit get a non-zero end position.
    n_reads <- length(subject_set)
    result <- integer(n_reads)

    n_matches <- lengths(matches)
    single_idx <- which(n_matches == 1L)

    if (length(single_idx) > 0L) {
        # Unlist just the single-match entries and pull their ends in one shot.
        result[single_idx] <- end(unlist(matches[single_idx], use.names = FALSE))
    }

    return(result)
}


#' Extract the cell barcode and UMI that immediately follow the adapter
#'
#' @param reads ShortReadQ object.
#' @param read1_end Integer vector with the adapter end coordinate per read.
#' @param len_CBC,len_UMI Lengths of the cell barcode and UMI.
#' @return List with character vectors \code{CBC} and \code{UMI}.
#' @noRd
.extract_cbc_umi <- function(reads, read1_end, len_CBC, len_UMI) {
    seqs <- sread(reads)

    # Vectorized subseq on a DNAStringSet extracts every barcode/UMI at once:
    # the CBC sits right after the adapter and the UMI right after the CBC.
    cbc_views <- subseq(seqs, start = read1_end + 1L, width = len_CBC)
    umi_views <- subseq(seqs, start = read1_end + len_CBC + 1L, width = len_UMI)

    list(
        CBC = as.character(cbc_views),
        UMI = as.character(umi_views)
    )
}


#' Call the base at an ROI using a fast two-pass strategy
#'
#' Pass 1 handles the common case where a read matches the ROI perfectly except
#' for the unknown base at the \code{N} position, extracting that base directly.
#' Pass 2 only processes the leftover reads (extra mismatches or indels) with an
#' expensive local pairwise alignment, memoised through a shared cache.
#'
#' @param reads ShortReadQ object for the current chunk.
#' @param roi_seq DNAString ROI containing a single \code{N} at the SNV position.
#' @param roi_mismatch Maximum mismatches tolerated in pass 2.
#' @param n_cores Cores used for the parallel alignment step.
#' @param cache Environment used to memoise alignment results across chunks.
#' @param n_pos_in_roi Unused; kept for backward compatibility of the signature.
#' @return Character vector (one per read) with the base at the ROI, or NA.
#' @noRd
.detect_roi_snv_optimized <- function(reads, roi_seq, roi_mismatch, n_cores,
                                      cache = NULL, n_pos_in_roi = NULL) {
    seqs <- sread(reads)
    n_reads <- length(seqs)

    # Pre-allocate the result as all-NA and fill in as calls are made.
    snv_calls <- character(n_reads)
    snv_calls[] <- NA_character_

    # Position of the unknown base (N) within the ROI, computed once.
    roi_seq_str <- as.character(roi_seq)
    n_pos_roi <- regexpr("N", roi_seq_str, fixed = TRUE)[1]
    roi_len <- nchar(roi_seq_str)

    # -------------------------------------------------------------------------
    # PASS 1: exact matching (only the N position may differ)
    # These reads match the ROI perfectly apart from the unknown base at N.
    # -------------------------------------------------------------------------
    matches_exact <- vmatchPattern2(roi_seq, seqs,
        max.mismatch = 1L,
        min.mismatch = 0L,
        fixed = TRUE,
        with.indels = FALSE # exact match: no indels
    )

    n_matches_exact <- lengths(matches_exact)
    exact_match_idx <- which(n_matches_exact == 1L)

    if (length(exact_match_idx) > 0L) {
        # Pull the matched region straight out of each read as a DNAStringSet.
        exact_seqs <- seqs[exact_match_idx]
        exact_ranges <- matches_exact[exact_match_idx]

        # extractAt returns a DNAStringSetList; unlist to a flat DNAStringSet.
        matched_seqs <- unlist(extractAt(exact_seqs, at = exact_ranges))

        # The matched region has the same layout as the ROI, so the called base
        # is simply the character at the N position.
        bases_at_n <- substr(as.character(matched_seqs), n_pos_roi, n_pos_roi)
        snv_calls[exact_match_idx] <- bases_at_n
    }

    # -------------------------------------------------------------------------
    # PASS 2: reads that did not match exactly need alignment
    # -------------------------------------------------------------------------
    remaining_idx <- which(is.na(snv_calls))

    if (length(remaining_idx) == 0L) {
        return(snv_calls)
    }

    # Re-match the leftovers with more mismatch tolerance and indels enabled.
    matches_fuzzy <- vmatchPattern2(roi_seq, seqs[remaining_idx],
        max.mismatch = roi_mismatch + 1L,
        min.mismatch = 0L,
        fixed = TRUE,
        with.indels = TRUE
    )

    n_matches_fuzzy <- lengths(matches_fuzzy)
    single_match_local <- which(n_matches_fuzzy == 1L)

    if (length(single_match_local) == 0L) {
        return(snv_calls)
    }

    # Translate local (within-remaining) indices back to original read indices.
    single_match_idx <- remaining_idx[single_match_local]

    # The alignment step is the expensive part, so run it in parallel.
    results <- mclapply(seq_along(single_match_local), function(j) {
        local_idx <- single_match_local[j]
        orig_idx <- single_match_idx[j]

        m <- matches_fuzzy[[local_idx]]
        pat_start <- start(m)
        pat_end <- end(m)

        # Extract the matched window from the original read.
        pat_seq <- as.character(subseq(seqs[[orig_idx]], start = pat_start, end = pat_end))

        # Reuse a previously computed call for an identical window if we have it.
        cached_result <- get0(pat_seq, envir = cache, inherits = FALSE)
        if (!is.null(cached_result)) {
            return(cached_result)
        }

        # Local alignment of the window against the ROI to place the N position.
        alignment <- tryCatch(
            {
                pwalign::pairwiseAlignment(DNAString(pat_seq), roi_seq,
                    type = "local",
                    gapOpening = 5,
                    gapExtension = 2
                )
            },
            error = function(e) NULL
        )

        if (is.null(alignment)) {
            return(NA_character_)
        }

        # Aligned strings let us find where the ROI's N landed in the read.
        aligned_pattern_str <- as.character(pwalign::pattern(alignment))
        aligned_reference_str <- as.character(pwalign::subject(alignment))

        n_pos <- regexpr("N", aligned_reference_str, fixed = TRUE)[1]

        if (n_pos > 0 && n_pos <= nchar(aligned_pattern_str)) {
            base_at_N <- substr(aligned_pattern_str, n_pos, n_pos)
            # Memoise so identical windows in later chunks are free.
            if (!is.null(cache)) {
                cache[[pat_seq]] <- base_at_N
            }
            return(base_at_N)
        } else {
            return(NA_character_)
        }
    }, mc.cores = n_cores)

    # Write the aligned calls back into the full-length result vector.
    snv_calls[single_match_idx] <- unlist(results)

    return(snv_calls)
}


#' Detect a poly(N) run (e.g. poly-T tail) in a set of sequences
#'
#' @param seqs DNAStringSet of the search windows.
#' @param polyN_base Nucleotide expected to repeat (e.g. "T").
#' @param polyN_length Required run length.
#' @param polyN_mismatch Mismatches tolerated within the run.
#' @return Logical vector: TRUE where the base count meets the threshold.
#' @noRd
.count_polyN_optimized <- function(seqs, polyN_base, polyN_length, polyN_mismatch) {
    min_required <- polyN_length - polyN_mismatch

    # alphabetFrequency is highly optimised; count the target base per window and
    # accept windows where it occurs often enough to be a poly(N) tail.
    base_matrix <- alphabetFrequency(seqs, baseOnly = TRUE)
    base_counts <- base_matrix[, polyN_base]

    return(base_counts >= min_required)
}


#' Process one FASTQ chunk into a table of SNV calls
#'
#' Applies the full per-read pipeline (length filter, Read 1 orientation and
#' chimera filter, poly(T) filter, CBC/UMI extraction, ROI base calling) and
#' returns a tidy data.table for the reads that survive.
#'
#' @return A \code{data.table} of surviving reads, or an empty one.
#' @noRd
.process_chunk <- function(chunk_files, read1_seq, read1_mismatch, read1_within_range,
                           polyN_base, polyN_length, polyN_mismatch,
                           major_roi_seq, major_roi_mismatch,
                           major_roi_mut, major_roi_wt,
                           minor_roi_seq, minor_roi_mismatch,
                           minor_roi_mut, minor_roi_wt,
                           len_CBC, len_UMI,
                           n_cores, major_snv_cache, minor_snv_cache,
                           use_minor_roi = TRUE, use_major_labeling = FALSE,
                           use_minor_labeling = FALSE) {
    n_input <- length(chunk_files)

    # Drop reads too short to possibly contain adapter + CBC + UMI + tail.
    min_length <- read1_within_range + len_CBC + len_UMI + polyN_length
    valid_idx <- width(chunk_files) >= min_length
    reads <- chunk_files[valid_idx]
    n_length_filter <- length(reads)

    if (n_length_filter == 0L) {
        message("  Chunk: No reads passed length filter")
        return(data.table())
    }

    # Grab the sequences once and reuse them throughout the chunk.
    seqs <- sread(reads)

    # Only the read start can contain the adapter, so search a trimmed window.
    trimmed_fwd <- subseq(seqs, start = 1L, end = read1_within_range)

    # Minimal per-read table carrying the index, id and full read width.
    dt <- data.table(
        index = seq_len(n_length_filter),
        id = as.character(id(reads)),
        width = width(reads)
    )

    # Forward-orientation adapter search.
    dt[, read1_fwd_end := .sequence_detection_vectorized(read1_seq, trimmed_fwd, read1_mismatch)]
    rm(trimmed_fwd) # free memory early

    # Reverse-complement search catches reads sequenced in the other orientation.
    seqs_rc <- reverseComplement(seqs)
    trimmed_rev <- subseq(seqs_rc, start = 1L, end = read1_within_range)

    dt[, read1_rev_end := .sequence_detection_vectorized(read1_seq, trimmed_rev, read1_mismatch)]
    rm(trimmed_rev)

    # A valid read has the adapter in exactly one orientation (XOR). Finding it
    # in both is a chimera and is discarded.
    dt[, `:=`(
        is_valid_read1 = xor(read1_fwd_end > 0L, read1_rev_end > 0L),
        is_dual_read = (read1_fwd_end > 0L) & (read1_rev_end > 0L)
    )]

    n_chimeric <- sum(dt$is_dual_read)
    n_valid_read1 <- sum(dt$is_valid_read1)

    if (n_valid_read1 == 0L) {
        message(sprintf("  Chunk: %d input -> %d length OK -> 0 valid Read1", n_input, n_length_filter))
        return(data.table())
    }

    # Split valid reads by the orientation in which the adapter was found.
    fwd_idx <- dt[is_valid_read1 == TRUE & read1_fwd_end > 0L, index]
    rev_idx <- dt[is_valid_read1 == TRUE & read1_rev_end > 0L, index]

    # Assemble one ShortReadQ of correctly-oriented reads. Reverse reads are
    # rebuilt from their reverse complement so downstream coordinates are uniform.
    if (length(rev_idx) > 0L) {
        # append() is required for ShortReadQ (c() would return a plain list).
        reads_fwd <- reads[fwd_idx]
        reads_rev <- ShortReadQ(
            sread = seqs_rc[rev_idx],
            quality = quality(reads[rev_idx]),
            id = id(reads[rev_idx])
        )
        reads_valid <- BiocGenerics::append(reads_fwd, reads_rev)
    } else {
        reads_valid <- reads[fwd_idx]
    }
    rm(seqs, seqs_rc, reads) # free memory
    gc(verbose = FALSE)

    # Matching table for the valid reads, in the same fwd-then-rev row order.
    dt_valid <- rbind(
        dt[is_valid_read1 == TRUE & read1_fwd_end > 0L],
        dt[is_valid_read1 == TRUE & read1_rev_end > 0L]
    )
    # For forward reads only the fwd end is non-zero (and vice versa), so summing
    # collapses the two columns into a single adapter end coordinate.
    dt_valid[, `:=`(
        read1_end_combined = read1_fwd_end + read1_rev_end,
        strand = fifelse(read1_fwd_end > 0L, "+", "-")
    )]
    rm(dt)

    # Look for the poly(T) tail in a window just after the UMI.
    polyT_search_window <- 50L
    polyT_search_start <- dt_valid$read1_end_combined + len_CBC + len_UMI + 1L
    polyT_search_end <- pmin(polyT_search_start + polyT_search_window, dt_valid$width)

    seqs_valid <- sread(reads_valid)
    reads_polyN_search <- subseq(seqs_valid, start = polyT_search_start, end = polyT_search_end)

    dt_valid[, polyN_exists := .count_polyN_optimized(reads_polyN_search, polyN_base, polyN_length, polyN_mismatch)]
    rm(reads_polyN_search)

    # Keep only reads with a detectable poly(T) tail.
    keep_idx <- which(dt_valid$polyN_exists)
    if (length(keep_idx) == 0L) {
        message(sprintf(
            "  Chunk: %d input -> %d length OK -> %d Read1 OK -> 0 polyT",
            n_input, n_length_filter, n_valid_read1
        ))
        return(data.table())
    }

    reads_final <- reads_valid[keep_idx]
    dt_final <- dt_valid[keep_idx]
    n_polyT_pass <- length(keep_idx)
    rm(reads_valid, dt_valid, seqs_valid)

    # Extract the cell barcode and UMI now that reads are oriented and filtered.
    barcodes <- .extract_cbc_umi(reads_final, dt_final$read1_end_combined, len_CBC, len_UMI)

    # Call the base at the major ROI (the mutation of interest).
    major_snv_base <- .detect_roi_snv_optimized(reads_final, major_roi_seq, major_roi_mismatch,
        n_cores,
        cache = major_snv_cache
    )

    # Call the base at the minor ROI (a nearby SNP used for phasing) if provided.
    if (use_minor_roi) {
        minor_snv_base <- .detect_roi_snv_optimized(reads_final, minor_roi_seq, minor_roi_mismatch,
            n_cores,
            cache = minor_snv_cache
        )
    } else {
        minor_snv_base <- rep(NA_character_, length(reads_final))
    }

    # Attach barcodes and called bases to the table by reference.
    dt_final[, `:=`(
        CBC = barcodes$CBC,
        UMI = barcodes$UMI,
        major_base = major_snv_base,
        minor_base = minor_snv_base
    )]

    # Optionally translate the called base into a MUT/WT label for the major ROI.
    if (use_major_labeling) {
        dt_final[, major_status := fcase(
            major_base == major_roi_mut, "MUT",
            major_base == major_roi_wt, "WT",
            default = NA_character_
        )]
    } else {
        dt_final[, major_status := NA_character_]
    }

    # Same optional MUT/WT labelling for the minor ROI.
    if (use_minor_roi && use_minor_labeling) {
        dt_final[, minor_status := fcase(
            minor_base == minor_roi_mut, "MUT",
            minor_base == minor_roi_wt, "WT",
            default = NA_character_
        )]
    } else {
        dt_final[, minor_status := NA_character_]
    }

    n_major_detected <- sum(!is.na(major_snv_base))

    message(sprintf(
        "  Chunk: %d -> %d len -> %d R1 -> %d chimeric -> %d polyT -> %d SNV (%.0f%%)",
        n_input, n_length_filter, n_valid_read1, n_chimeric,
        n_polyT_pass, n_major_detected, 100 * n_major_detected / n_polyT_pass
    ))

    # Present the most useful columns first and rename id -> read_id.
    setcolorder(dt_final, c(
        "id", "CBC", "UMI", "strand",
        "major_base", "major_status",
        "minor_base", "minor_status"
    ))

    setnames(dt_final, "id", "read_id")

    return(dt_final)
}


# --- Main function -----------------------------------------------------------

#' Detect SNVs in single-cell long-read sequencing reads
#'
#' Streams a FASTQ file in chunks and, for every read, detects the Read 1
#' adapter, extracts the cell barcode and UMI, requires a poly(T) tail, and calls
#' the base at a major region of interest (the mutation) and an optional minor
#' region of interest (a SNP used later for phasing). Results are written to a
#' CSV that feeds \code{\link{summarize_snv}}.
#'
#' @param session_name Prefix used for output files.
#' @param path_input_file Path to a single FASTQ file (\code{.fastq} or
#'   \code{.fastq.gz}).
#' @param path_output_folder Existing directory to write output files into.
#' @param read1_seq Adapter/Read 1 sequence that precedes the cell barcode.
#' @param read1_mismatch Allowed mismatches when matching \code{read1_seq}.
#' @param read1_within_range Size of the read-start window searched for the
#'   adapter (default: 70).
#' @param polyN_base Nucleotide of the tail to require (default: "T").
#' @param polyN_length Required poly(N) run length (default: 10).
#' @param polyN_mismatch Mismatches tolerated within the poly(N) run (default: 2).
#' @param major_roi_seq ROI sequence with a single \code{N} at the major SNV
#'   position.
#' @param major_roi_mismatch Allowed mismatches when matching \code{major_roi_seq}.
#' @param major_roi_mut Base that denotes the mutation in the major ROI
#'   (optional; enables MUT/WT labelling).
#' @param major_roi_wt Base that denotes wild-type in the major ROI (optional).
#' @param minor_roi_seq ROI sequence with a single \code{N} at a SNP position
#'   used for phasing (optional).
#' @param minor_roi_mismatch Allowed mismatches for \code{minor_roi_seq} (optional).
#' @param minor_roi_mut Base that denotes the mutation in the minor ROI (optional).
#' @param minor_roi_wt Base that denotes wild-type in the minor ROI (optional).
#' @param len_CBC Cell barcode length (default: 16).
#' @param len_UMI UMI length (default: 12).
#' @param n_cores Number of CPU cores for the parallel alignment step.
#' @param n_reads_per_chunk Reads processed per streamed chunk (default: 1e5).
#'
#' @return Invisibly, nothing useful is returned; the SNV table is written to
#'   \code{<session_name>_snv_table.csv} in \code{path_output_folder}.
#'
#' @seealso \code{\link{summarize_snv}}, \code{\link{flag_snv}}
#' @export
detect_snv <- function(session_name,
                       path_input_file,
                       path_output_folder,
                       read1_seq,
                       read1_mismatch,
                       read1_within_range = 70,
                       polyN_base = "T",
                       polyN_length = 10,
                       polyN_mismatch = 2,
                       major_roi_seq,
                       major_roi_mismatch,
                       major_roi_mut = NULL,
                       major_roi_wt = NULL,
                       minor_roi_seq = NULL,
                       minor_roi_mismatch = NULL,
                       minor_roi_mut = NULL,
                       minor_roi_wt = NULL,
                       len_CBC = 16,
                       len_UMI = 12,
                       n_cores,
                       n_reads_per_chunk = 1e5) {
    # --- Input validation ---
    stopifnot("`session_name` must be character" = is.character(session_name))
    stopifnot("`session_name` must not contain special characters" = grepl("^[a-zA-Z0-9_-]+$", session_name))
    stopifnot("`path_input_file` must be character" = is.character(path_input_file))
    stopifnot("`path_input_file` - file not found" = file.exists(path_input_file))
    stopifnot("`path_input_file` must be .fastq or .fastq.gz" = grepl("\\.fastq(\\.gz)?$", path_input_file))
    stopifnot("`path_output_folder` must be character" = is.character(path_output_folder))
    stopifnot("`path_output_folder` - path not found" = dir.exists(path_output_folder))
    stopifnot("`read1_seq` must contain only A,C,G,T" = grepl("^[GATC]+$", read1_seq))
    stopifnot("`read1_mismatch` must be integer" = is.numeric(read1_mismatch) && read1_mismatch %% 1 == 0)
    stopifnot("`major_roi_seq` must have exactly one N" = grepl("^[GATC]*N[GATC]*$", major_roi_seq))
    stopifnot("`major_roi_mismatch` must be integer" = is.numeric(major_roi_mismatch) && major_roi_mismatch %% 1 == 0)
    stopifnot("`n_cores` must be integer" = is.numeric(n_cores) && n_cores %% 1 == 0)

    # Major-ROI labelling is enabled only when both mut and wt bases are supplied.
    use_major_labeling <- !is.null(major_roi_mut) || !is.null(major_roi_wt)
    if (use_major_labeling) {
        stopifnot(
            "Both major_roi_mut and major_roi_wt must be provided together" =
                !is.null(major_roi_mut) && !is.null(major_roi_wt)
        )
        stopifnot("`major_roi_mut` must be single nucleotide" = grepl("^[ACGT]$", major_roi_mut))
        stopifnot("`major_roi_wt` must be single nucleotide" = grepl("^[ACGT]$", major_roi_wt))
        stopifnot("`major_roi_mut` and `major_roi_wt` must differ" = major_roi_mut != major_roi_wt)
        message("  Major ROI labeling: ", major_roi_mut, " = MUT, ", major_roi_wt, " = WT")
    }

    # The minor ROI (and its optional labelling) is entirely opt-in.
    use_minor_roi <- !is.null(minor_roi_seq)
    use_minor_labeling <- FALSE

    if (use_minor_roi) {
        stopifnot("`minor_roi_seq` must have exactly one N" = grepl("^[GATC]*N[GATC]*$", minor_roi_seq))
        stopifnot("`minor_roi_mismatch` required" = !is.null(minor_roi_mismatch))

        use_minor_labeling <- !is.null(minor_roi_mut) || !is.null(minor_roi_wt)
        if (use_minor_labeling) {
            stopifnot(
                "Both minor_roi_mut and minor_roi_wt must be provided together" =
                    !is.null(minor_roi_mut) && !is.null(minor_roi_wt)
            )
            message("  Minor ROI labeling: ", minor_roi_mut, " = MUT, ", minor_roi_wt, " = WT")
        }
    } else {
        message("  Minor ROI: not provided, skipping")
    }

    # Never request more than (available - 1) cores.
    available_cores <- parallel::detectCores()
    if (n_cores > available_cores - 1) {
        n_cores <- max(1L, available_cores - 1L)
        message("  n_cores adjusted to ", n_cores)
    }

    # Convert search sequences to DNAString once up front.
    read1_seq <- DNAString(read1_seq)
    major_roi_seq <- DNAString(major_roi_seq)

    if (use_minor_roi) {
        minor_roi_seq <- DNAString(minor_roi_seq)
    }

    # Alignment caches shared across chunks so repeated windows align only once.
    major_snv_cache <- new.env(hash = TRUE, parent = emptyenv())
    minor_snv_cache <- if (use_minor_roi) new.env(hash = TRUE, parent = emptyenv()) else NULL

    # Output CSV (appended to, chunk by chunk).
    out_file <- file.path(path_output_folder, paste0(session_name, "_snv_table.csv"))

    # Stream the FASTQ in chunks to keep memory bounded on large files.
    message(Sys.time(), " - Starting SNV detection")
    message("  Input: ", path_input_file)
    message("  Chunk size: ", format(n_reads_per_chunk, scientific = FALSE))

    fq_stream <- FastqStreamer(path_input_file, n = n_reads_per_chunk)
    on.exit(close(fq_stream), add = TRUE)

    chunk_id <- 1L
    total_reads <- 0L
    total_detected <- 0L

    repeat {
        fq_chunk <- yield(fq_stream)
        if (length(fq_chunk) == 0L) {
            message("All chunks processed.")
            break
        }

        results <- .process_chunk(
            fq_chunk,
            read1_seq, read1_mismatch, read1_within_range,
            polyN_base, polyN_length, polyN_mismatch,
            major_roi_seq, major_roi_mismatch,
            major_roi_mut, major_roi_wt,
            minor_roi_seq, minor_roi_mismatch,
            minor_roi_mut, minor_roi_wt,
            len_CBC, len_UMI, n_cores, major_snv_cache, minor_snv_cache,
            use_minor_roi, use_major_labeling, use_minor_labeling
        )

        # Append surviving reads; only the first chunk writes the header.
        if (nrow(results) > 0) {
            fwrite(results, out_file, append = (chunk_id > 1L))
            total_detected <- total_detected + nrow(results)
        }

        total_reads <- total_reads + length(fq_chunk)

        message(sprintf(
            "%s - Chunk %d: %d reads processed, %d total detected",
            Sys.time(), chunk_id, length(fq_chunk), total_detected
        ))

        chunk_id <- chunk_id + 1L
        gc(verbose = FALSE)
    }

    message("\n--- Complete ---")
    message("Total reads: ", format(total_reads, big.mark = ","))
    message("SNVs detected: ", format(total_detected, big.mark = ","))
    message("Cache hits (major): ", length(major_snv_cache))
    if (!is.null(minor_snv_cache)) {
        message("Cache hits (minor): ", length(minor_snv_cache))
    }
    message("Output: ", out_file)
}
