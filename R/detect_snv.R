# =============================================================================
# detect_snv.R - SNV Detection from Nanopore Sequencing Reads
# =============================================================================

library(ShortRead)
library(Biostrings)
library(BiocGenerics)
library(pwalign)
library(data.table)

# --- Internal helpers ---

# Vectorized pattern detection returning end position of single matches
.sequence_detection_vectorized <- function(pattern, subject_set, max_mismatch) {
  matches <- vmatchPattern2(pattern, subject_set,
                            max.mismatch = max_mismatch,
                            min.mismatch = 0,
                            fixed = TRUE,
                            with.indels = TRUE
  )
  n_reads <- length(subject_set)
  result <- integer(n_reads)
  n_matches <- lengths(matches)
  single_idx <- which(n_matches == 1L)
  if (length(single_idx) > 0L) {
    result[single_idx] <- end(unlist(matches[single_idx], use.names = FALSE))
  }
  return(result)
}


# Detailed pattern detection returning end position AND match width
# Used by enriched mode to classify Read1 error types (indels change width)
.sequence_detection_detailed <- function(pattern, subject_set, max_mismatch) {
  matches <- vmatchPattern2(pattern, subject_set,
                            max.mismatch = max_mismatch, min.mismatch = 0,
                            fixed = TRUE, with.indels = TRUE
  )
  n_reads <- length(subject_set)
  end_pos <- integer(n_reads)
  match_width <- integer(n_reads)
  matched_seq <- character(n_reads)
  matched_seq[] <- NA_character_
  n_matches <- lengths(matches)
  single_idx <- which(n_matches == 1L)
  if (length(single_idx) > 0L) {
    unlisted <- unlist(matches[single_idx], use.names = FALSE)
    end_pos[single_idx] <- end(unlisted)
    match_width[single_idx] <- width(unlisted)
    extracted <- Biostrings::extractAt(subject_set[single_idx], at = matches[single_idx])
    matched_seq[single_idx] <- as.character(unlist(extracted))
  }
  list(end_pos = end_pos, match_width = match_width, matched_seq = matched_seq)
}


# Extracts cell barcode (CBC) and UMI using vectorized XStringSet operations
.extract_cbc_umi <- function(reads, read1_end, len_CBC, len_UMI) {
  seqs <- ShortRead::sread(reads)
  cbc_views <- Biostrings::subseq(seqs, start = read1_end + 1L, width = len_CBC)
  umi_views <- Biostrings::subseq(seqs, start = read1_end + len_CBC + 1L, width = len_UMI)
  list(
    CBC = as.character(cbc_views),
    UMI = as.character(umi_views)
  )
}


# Two-pass ROI SNV detection
# Pass 1: exact IUPAC match (fixed=FALSE, 0 mismatches) - direct extraction
# Pass 2: fuzzy IUPAC match + pairwise alignment for remaining reads
# When enriched=TRUE, returns list with base, method, matched_seq, error counts, positions
.detect_roi_snv_optimized <- function(reads, roi_seq, roi_mismatch, n_cores,
                                      cache = NULL, n_pos_in_roi = NULL,
                                      enriched = FALSE) {
  seqs <- ShortRead::sread(reads)
  n_reads <- length(seqs)
  snv_calls <- character(n_reads)
  snv_calls[] <- NA_character_
  
  roi_seq_str <- as.character(roi_seq)
  roi_len <- nchar(roi_seq_str)
  n_pos_roi <- regexpr("N", roi_seq_str, fixed = TRUE)[1]
  
  # Enriched output vectors
  if (enriched) {
    roi_method <- character(n_reads)
    roi_method[] <- NA_character_
    roi_matched_seq <- character(n_reads)
    roi_matched_seq[] <- NA_character_
    roi_n_sub <- integer(n_reads)
    roi_n_sub[] <- NA_integer_
    roi_n_del <- integer(n_reads)
    roi_n_del[] <- NA_integer_
    roi_n_ins <- integer(n_reads)
    roi_n_ins[] <- NA_integer_
    roi_start_pos <- integer(n_reads)
    roi_start_pos[] <- NA_integer_
    roi_end_pos <- integer(n_reads)
    roi_end_pos[] <- NA_integer_
    roi_n_offset <- integer(n_reads)
    roi_n_offset[] <- NA_integer_
  }
  
  # Pass 1: exact IUPAC match (N and other degenerate codes are free)
  matches_exact <- vmatchPattern2(roi_seq, seqs,
                                  max.mismatch = 0L, min.mismatch = 0L,
                                  fixed = FALSE, with.indels = FALSE
  )
  n_matches_exact <- lengths(matches_exact)
  exact_match_idx <- which(n_matches_exact == 1L)
  
  if (length(exact_match_idx) > 0L) {
    exact_seqs <- seqs[exact_match_idx]
    exact_ranges <- matches_exact[exact_match_idx]
    matched_seqs <- unlist(Biostrings::extractAt(exact_seqs, at = exact_ranges))
    matched_seqs_char <- as.character(matched_seqs)
    bases_at_n <- substr(matched_seqs_char, n_pos_roi, n_pos_roi)
    snv_calls[exact_match_idx] <- bases_at_n
    
    if (enriched) {
      roi_method[exact_match_idx] <- "exact"
      roi_matched_seq[exact_match_idx] <- matched_seqs_char
      # Pass 1 = 0 mismatches with IUPAC expansion, so 0 errors by definition
      roi_n_sub[exact_match_idx] <- 0L
      roi_n_del[exact_match_idx] <- 0L
      roi_n_ins[exact_match_idx] <- 0L
      # ROI match positions (no indels, so N offset is fixed)
      exact_ranges_flat <- unlist(exact_ranges, use.names = FALSE)
      roi_start_pos[exact_match_idx] <- start(exact_ranges_flat)
      roi_end_pos[exact_match_idx] <- end(exact_ranges_flat)
      roi_n_offset[exact_match_idx] <- n_pos_roi
    }
  }
  
  # Pass 2: fuzzy match + alignment
  remaining_idx <- which(is.na(snv_calls))
  if (length(remaining_idx) == 0L) {
    if (enriched) {
      return(list(base = snv_calls, method = roi_method,
                  matched_seq = roi_matched_seq,
                  n_substitutions = roi_n_sub, n_deletions = roi_n_del, n_insertions = roi_n_ins,
                  roi_start = roi_start_pos, roi_end = roi_end_pos, n_offset = roi_n_offset))
    }
    return(snv_calls)
  }
  
  matches_fuzzy <- vmatchPattern2(roi_seq, seqs[remaining_idx],
                                  max.mismatch = roi_mismatch, min.mismatch = 0L,
                                  fixed = FALSE, with.indels = TRUE
  )
  n_matches_fuzzy <- lengths(matches_fuzzy)
  n_zero <- sum(n_matches_fuzzy == 0L)
  n_single <- sum(n_matches_fuzzy == 1L)
  n_multi <- sum(n_matches_fuzzy > 1L)
  message(sprintf("    Pass 2: %d remaining -> %d no match, %d single, %d multi-match",
                  length(remaining_idx), n_zero, n_single, n_multi))
  single_match_local <- which(n_matches_fuzzy == 1L)
  if (length(single_match_local) == 0L) {
    if (enriched) {
      return(list(base = snv_calls, method = roi_method,
                  matched_seq = roi_matched_seq,
                  n_substitutions = roi_n_sub, n_deletions = roi_n_del, n_insertions = roi_n_ins,
                  roi_start = roi_start_pos, roi_end = roi_end_pos, n_offset = roi_n_offset))
    }
    return(snv_calls)
  }
  
  single_match_idx <- remaining_idx[single_match_local]
  
  results <- parallel::mclapply(seq_along(single_match_local), function(j) {
    local_idx <- single_match_local[j]
    orig_idx <- single_match_idx[j]
    m <- matches_fuzzy[[local_idx]]
    match_w <- width(m)
    pat_seq <- as.character(Biostrings::subseq(seqs[[orig_idx]],
                                               start = start(m), end = end(m)))
    
    cached_result <- get0(pat_seq, envir = cache, inherits = FALSE)
    if (!is.null(cached_result)) return(cached_result)
    
    tryCatch({
      # Single padded global alignment — used for both base calling and error counting
      pad <- "ACGTACGTAC"
      padded_pat <- paste0(pad, pat_seq, pad)
      padded_ref <- paste0(pad, roi_seq_str, pad)
      aln <- pwalign::pairwiseAlignment(
        Biostrings::DNAString(padded_pat), Biostrings::DNAString(padded_ref),
        type = "global", gapOpening = 5, gapExtension = 2)
      
      # Find N position in aligned reference (within padded alignment)
      aligned_ref_str <- as.character(pwalign::alignedSubject(aln))
      aligned_pat_str <- as.character(pwalign::alignedPattern(aln))
      n_pos <- regexpr("N", aligned_ref_str, fixed = TRUE)[1]
      
      if (n_pos > 0 && n_pos <= nchar(aligned_pat_str)) {
        base_at_N <- substr(aligned_pat_str, n_pos, n_pos)
        if (base_at_N == "-") base_at_N <- NA_character_
        
        if (enriched) {
          n_sub <- pwalign::nmismatch(aln)
          ins <- pwalign::insertion(aln)[[1]]
          del <- pwalign::deletion(aln)[[1]]
          
          # N offset: count non-gap chars in aligned pattern up to n_pos
          if (!is.na(base_at_N)) {
            aln_pat_chars <- strsplit(aligned_pat_str, "")[[1]]
            n_non_gap <- sum(aln_pat_chars[seq_len(n_pos)] != "-")
            # Subtract pad length to get offset within original pat_seq
            n_offset_val <- as.integer(n_non_gap - nchar(pad))
          } else {
            n_offset_val <- NA_integer_
          }
          
          result <- list(base = base_at_N, matched_seq = pat_seq,
                         n_substitutions = n_sub,
                         n_deletions = sum(width(del)),
                         n_insertions = sum(width(ins)),
                         n_offset_in_match = n_offset_val)
          return(result)
        }
        
        return(base_at_N)
      }
      if (enriched) return(list(base = NA_character_, matched_seq = pat_seq,
                                n_substitutions = NA_integer_, n_deletions = NA_integer_, n_insertions = NA_integer_,
                                n_offset_in_match = NA_integer_))
      NA_character_
    }, error = function(e) {
      if (enriched) return(list(base = NA_character_, matched_seq = pat_seq,
                                n_substitutions = NA_integer_, n_deletions = NA_integer_, n_insertions = NA_integer_,
                                n_offset_in_match = NA_integer_))
      NA_character_
    })
  }, mc.cores = n_cores)
  
  # Process results and populate cache in one pass
  if (enriched) {
    fuzzy_ranges_flat <- unlist(matches_fuzzy[single_match_local], use.names = FALSE)
    roi_start_pos[single_match_idx] <- start(fuzzy_ranges_flat)
    roi_end_pos[single_match_idx] <- end(fuzzy_ranges_flat)
  }
  
  for (k in seq_along(single_match_local)) {
    r <- results[[k]]
    idx <- single_match_idx[k]
    
    if (enriched && is.list(r)) {
      snv_calls[idx] <- r$base
      roi_method[idx] <- "alignment"
      roi_matched_seq[idx] <- r$matched_seq
      roi_n_sub[idx] <- r$n_substitutions
      roi_n_del[idx] <- r$n_deletions
      roi_n_ins[idx] <- r$n_insertions
      roi_n_offset[idx] <- r$n_offset_in_match
      # Populate cache using matched_seq already in the result
      if (!is.null(cache) && !is.null(r$matched_seq)) {
        if (!exists(r$matched_seq, envir = cache, inherits = FALSE)) {
          cache[[r$matched_seq]] <- r
        }
      }
    } else if (!enriched) {
      # roi_method only exists in enriched mode; just record the base
      if (is.character(r)) {
        snv_calls[idx] <- r
      }
      # Cache for non-enriched: extract pat_seq
      if (!is.null(cache)) {
        local_idx <- single_match_local[k]
        orig_idx <- single_match_idx[k]
        m <- matches_fuzzy[[local_idx]]
        ps <- as.character(Biostrings::subseq(seqs[[orig_idx]],
                                              start = start(m), end = end(m)))
        if (!exists(ps, envir = cache, inherits = FALSE)) {
          cache[[ps]] <- r
        }
      }
    } else {
      # enriched but r is not a list (error case)
      snv_calls[idx] <- if (is.character(r)) r else NA_character_
      roi_method[idx] <- "alignment"
    }
  }
  
  if (enriched) {
    return(list(base = snv_calls, method = roi_method,
                matched_seq = roi_matched_seq,
                n_substitutions = roi_n_sub, n_deletions = roi_n_del, n_insertions = roi_n_ins,
                roi_start = roi_start_pos, roi_end = roi_end_pos, n_offset = roi_n_offset))
  }
  
  # snv_calls was already populated in the loop above
  return(snv_calls)
}


# Optimized polyN detection using alphabetFrequency
.count_polyN_optimized <- function(seqs, polyN_base, polyN_length, polyN_mismatch) {
  min_required <- polyN_length - polyN_mismatch
  base_matrix <- Biostrings::alphabetFrequency(seqs, baseOnly = TRUE)
  base_counts <- base_matrix[, polyN_base]
  return(base_counts >= min_required)
}


# Counts substitutions, deletions, and insertions between two sequences
# Uses padded global alignment to correctly detect terminal errors
# Returns error counts via pwalign built-in functions (nmismatch, insertion, deletion)
.count_errors_padded <- function(read_seq, ref_seq) {
  pad <- "ACGTACGTAC"
  padded_read <- paste0(pad, read_seq, pad)
  padded_ref <- paste0(pad, ref_seq, pad)
  aln <- pwalign::pairwiseAlignment(
    Biostrings::DNAString(padded_read), Biostrings::DNAString(padded_ref),
    type = "global", gapOpening = 5, gapExtension = 2)
  n_sub <- pwalign::nmismatch(aln)
  ins <- pwalign::insertion(aln)[[1]]
  del <- pwalign::deletion(aln)[[1]]
  list(n_substitutions = n_sub, n_deletions = sum(width(del)), n_insertions = sum(width(ins)))
}


# Processes a single chunk of FASTQ reads
.process_chunk <- function(chunk_files, read1_seq, read1_mismatch, read1_within_range,
                           polyN_base, polyN_length, polyN_mismatch,
                           major_roi_seq, major_roi_mismatch,
                           major_roi_mut, major_roi_wt,
                           minor_roi_seq, minor_roi_mismatch,
                           minor_roi_mut, minor_roi_wt,
                           len_CBC, len_UMI,
                           n_cores, major_snv_cache, minor_snv_cache,
                           read1_error_cache = NULL,
                           use_minor_roi = TRUE, use_major_labeling = FALSE,
                           use_minor_labeling = FALSE,
                           enriched = FALSE) {
  n_input <- length(chunk_files)
  min_length <- read1_within_range + len_CBC + len_UMI +
    (if (is.null(polyN_base)) 0L else polyN_length)
  valid_idx <- width(chunk_files) >= min_length
  reads <- chunk_files[valid_idx]
  n_length_filter <- length(reads)
  
  if (n_length_filter == 0L) {
    message("  Chunk: No reads passed length filter")
    return(list(data = data.table::data.table(), stats = list(
      n_input = n_input, n_length_filter = 0L, n_valid_read1 = 0L,
      n_chimeric = 0L, n_polyT_pass = 0L, n_major_detected = 0L, n_minor_detected = 0L)))
  }
  
  seqs <- ShortRead::sread(reads)
  trimmed_fwd <- Biostrings::subseq(seqs, start = 1L, end = read1_within_range)
  
  dt <- data.table::data.table(
    index = seq_len(n_length_filter),
    id = as.character(ShortRead::id(reads)),
    width = width(reads)
  )
  
  if (enriched) {
    fwd_detail <- .sequence_detection_detailed(read1_seq, trimmed_fwd, read1_mismatch)
    dt[, read1_fwd_end := fwd_detail$end_pos]
    dt[, read1_fwd_width := fwd_detail$match_width]
  } else {
    dt[, read1_fwd_end := .sequence_detection_vectorized(read1_seq, trimmed_fwd, read1_mismatch)]
  }
  rm(trimmed_fwd)
  
  seqs_rc <- Biostrings::reverseComplement(seqs)
  trimmed_rev <- Biostrings::subseq(seqs_rc, start = 1L, end = read1_within_range)
  if (enriched) {
    rev_detail <- .sequence_detection_detailed(read1_seq, trimmed_rev, read1_mismatch)
    dt[, read1_rev_end := rev_detail$end_pos]
    dt[, read1_rev_width := rev_detail$match_width]
  } else {
    dt[, read1_rev_end := .sequence_detection_vectorized(read1_seq, trimmed_rev, read1_mismatch)]
  }
  rm(trimmed_rev)
  
  dt[, `:=`(
    is_valid_read1 = xor(read1_fwd_end > 0L, read1_rev_end > 0L),
    is_dual_read = (read1_fwd_end > 0L) & (read1_rev_end > 0L)
  )]
  
  n_chimeric <- sum(dt$is_dual_read)
  n_valid_read1 <- sum(dt$is_valid_read1)
  
  if (n_valid_read1 == 0L) {
    message(sprintf("  Chunk: %d input -> %d length OK -> 0 valid Read1", n_input, n_length_filter))
    return(list(data = data.table::data.table(), stats = list(
      n_input = n_input, n_length_filter = n_length_filter, n_valid_read1 = 0L,
      n_chimeric = n_chimeric, n_polyT_pass = 0L, n_major_detected = 0L, n_minor_detected = 0L)))
  }
  
  fwd_idx <- dt[is_valid_read1 == TRUE & read1_fwd_end > 0L, index]
  rev_idx <- dt[is_valid_read1 == TRUE & read1_rev_end > 0L, index]
  
  if (length(rev_idx) > 0L) {
    reads_fwd <- reads[fwd_idx]
    reads_rev <- ShortRead::ShortReadQ(
      sread = seqs_rc[rev_idx],
      quality = Biostrings::quality(reads[rev_idx]),
      id = ShortRead::id(reads[rev_idx])
    )
    reads_valid <- BiocGenerics::append(reads_fwd, reads_rev)
  } else {
    reads_valid <- reads[fwd_idx]
  }
  rm(seqs, seqs_rc, reads)
  gc(verbose = FALSE)
  
  dt_valid <- rbind(
    dt[is_valid_read1 == TRUE & read1_fwd_end > 0L],
    dt[is_valid_read1 == TRUE & read1_rev_end > 0L]
  )
  dt_valid[, `:=`(
    read1_end_combined = read1_fwd_end + read1_rev_end,
    strand = data.table::fifelse(read1_fwd_end > 0L, "+", "-")
  )]
  rm(dt)
  
  seqs_valid <- ShortRead::sread(reads_valid)
  
  # Poly(N) tail filter. Skipped entirely when polyN_base is NULL
  # (e.g. 5' GEX libraries, where the TSO rather than a poly(T) follows the UMI).
  if (!is.null(polyN_base)) {
    polyT_search_window <- 50L
    polyT_search_start <- dt_valid$read1_end_combined + len_CBC + len_UMI + 1L
    polyT_search_end <- pmin(polyT_search_start + polyT_search_window, dt_valid$width)
    
    reads_polyN_search <- Biostrings::subseq(seqs_valid, start = polyT_search_start, end = polyT_search_end)
    dt_valid[, polyN_exists := .count_polyN_optimized(reads_polyN_search, polyN_base, polyN_length, polyN_mismatch)]
    rm(reads_polyN_search)
  } else {
    dt_valid[, polyN_exists := TRUE]
  }
  
  keep_idx <- which(dt_valid$polyN_exists)
  if (length(keep_idx) == 0L) {
    message(sprintf("  Chunk: %d input -> %d length OK -> %d Read1 OK -> 0 polyT",
                    n_input, n_length_filter, n_valid_read1))
    return(list(data = data.table::data.table(), stats = list(
      n_input = n_input, n_length_filter = n_length_filter, n_valid_read1 = n_valid_read1,
      n_chimeric = n_chimeric, n_polyT_pass = 0L, n_major_detected = 0L, n_minor_detected = 0L)))
  }
  
  reads_final <- reads_valid[keep_idx]
  dt_final <- dt_valid[keep_idx]
  n_polyT_pass <- length(keep_idx)
  rm(reads_valid, dt_valid, seqs_valid)
  
  barcodes <- .extract_cbc_umi(reads_final, dt_final$read1_end_combined, len_CBC, len_UMI)
  
  # ROI detection (enriched or standard)
  major_roi_result <- .detect_roi_snv_optimized(reads_final, major_roi_seq, major_roi_mismatch,
                                                n_cores, cache = major_snv_cache, enriched = enriched)
  
  if (use_minor_roi) {
    minor_roi_result <- .detect_roi_snv_optimized(reads_final, minor_roi_seq, minor_roi_mismatch,
                                                  n_cores, cache = minor_snv_cache, enriched = enriched)
  }
  
  if (enriched) {
    major_snv_base <- major_roi_result$base
    minor_snv_base <- if (use_minor_roi) minor_roi_result$base else rep(NA_character_, length(reads_final))
  } else {
    major_snv_base <- major_roi_result
    minor_snv_base <- if (use_minor_roi) minor_roi_result else rep(NA_character_, length(reads_final))
  }
  
  dt_final[, `:=`(
    CBC = barcodes$CBC,
    UMI = barcodes$UMI,
    major_base = major_snv_base,
    minor_base = minor_snv_base
  )]
  
  if (use_major_labeling) {
    dt_final[, major_status := data.table::fcase(
      major_base == major_roi_mut, "MUT",
      major_base == major_roi_wt, "WT",
      default = NA_character_
    )]
  } else {
    dt_final[, major_status := NA_character_]
  }
  
  if (use_minor_roi && use_minor_labeling) {
    dt_final[, minor_status := data.table::fcase(
      minor_base == minor_roi_mut, "MUT",
      minor_base == minor_roi_wt, "WT",
      default = NA_character_
    )]
  } else {
    dt_final[, minor_status := NA_character_]
  }
  
  # --- Enriched columns ---
  if (enriched) {
    read1_pattern_len <- nchar(as.character(read1_seq))
    
    # Read length
    dt_final[, read_length := width]
    
    # Read1 position (already in read1_end_combined)
    dt_final[, read1_pos := read1_end_combined]
    
    # Read1 match width and error counts
    r1_width <- data.table::fifelse(
      dt_final$strand == "+",
      dt_final$read1_fwd_width,
      dt_final$read1_rev_width
    )
    dt_final[, read1_match_width := r1_width]
    
    # Get matched sequences for the active strand
    r1_matched <- character(nrow(dt_final))
    is_fwd <- dt_final$strand == "+"
    r1_matched[is_fwd] <- fwd_detail$matched_seq[dt_final$index[is_fwd]]
    r1_matched[!is_fwd] <- rev_detail$matched_seq[dt_final$index[!is_fwd]]
    
    # Count errors via pairwise alignment — deduplicate + cache across chunks
    read1_pat_str <- as.character(read1_seq)
    r1_n_sub <- integer(nrow(dt_final))
    r1_n_del <- integer(nrow(dt_final))
    r1_n_ins <- integer(nrow(dt_final))
    match_idx <- which(!is.na(r1_matched))
    
    if (length(match_idx) > 0L) {
      # Deduplicate: only align unique sequences
      unique_seqs <- unique(r1_matched[match_idx])
      
      # Check cache for already-computed sequences
      cached_mask <- logical(length(unique_seqs))
      cached_results <- vector("list", length(unique_seqs))
      if (!is.null(read1_error_cache)) {
        for (u in seq_along(unique_seqs)) {
          cr <- get0(unique_seqs[u], envir = read1_error_cache, inherits = FALSE)
          if (!is.null(cr)) {
            cached_mask[u] <- TRUE
            cached_results[[u]] <- cr
          }
        }
      }
      
      # Align only uncached unique sequences
      to_align <- which(!cached_mask)
      if (length(to_align) > 0L) {
        aln_results <- parallel::mclapply(to_align, function(u) {
          tryCatch(
            .count_errors_padded(unique_seqs[u], read1_pat_str),
            error = function(e) list(n_substitutions = NA_integer_,
                                     n_deletions = NA_integer_, n_insertions = NA_integer_))
        }, mc.cores = n_cores)
        # Store results and populate cache in parent process
        for (k in seq_along(to_align)) {
          u <- to_align[k]
          r <- aln_results[[k]]
          if (!is.list(r)) r <- list(n_substitutions = NA_integer_,
                                     n_deletions = NA_integer_, n_insertions = NA_integer_)
          cached_results[[u]] <- r
          if (!is.null(read1_error_cache)) {
            read1_error_cache[[unique_seqs[u]]] <- r
          }
        }
      }
      
      # Build lookup vectors from unique results
      res_sub <- vapply(cached_results, function(r) r$n_substitutions, integer(1))
      res_del <- vapply(cached_results, function(r) r$n_deletions, integer(1))
      res_ins <- vapply(cached_results, function(r) r$n_insertions, integer(1))
      
      # Map results back to all reads via match (vectorized)
      lookup_idx <- match(r1_matched[match_idx], unique_seqs)
      r1_n_sub[match_idx] <- res_sub[lookup_idx]
      r1_n_del[match_idx] <- res_del[lookup_idx]
      r1_n_ins[match_idx] <- res_ins[lookup_idx]
    }
    
    no_match <- is.na(r1_matched)
    r1_n_sub[no_match] <- NA_integer_
    r1_n_del[no_match] <- NA_integer_
    r1_n_ins[no_match] <- NA_integer_
    
    dt_final[, `:=`(
      read1_n_substitutions = r1_n_sub,
      read1_n_deletions = r1_n_del,
      read1_n_insertions = r1_n_ins
    )]
    
    # Q-scores — vectorized using IntegerList operations
    qual_bstring <- Biostrings::quality(Biostrings::quality(reads_final))
    qual_int <- as(qual_bstring, "IntegerList") - 33L  # Phred+33 decoding
    qual_lens <- lengths(qual_int)
    
    # Mean Q-score: vectorized via sum/lengths
    dt_final[, mean_qscore := round(sum(qual_int) / qual_lens, 1)]
    
    # Read1 Q-score: extract sublist [1:read1_end] per read, then mean
    r1_ends <- pmin(dt_final$read1_end_combined, qual_lens)
    r1_ranges <- IRanges::IRanges(start = 1L, end = r1_ends)
    qual_r1 <- IRanges::extractList(unlist(qual_int, use.names = FALSE),
                                    IRanges::shift(r1_ranges, shift = c(0L, cumsum(qual_lens)[-length(qual_lens)])))
    dt_final[, read1_qscore := round(sum(qual_r1) / lengths(qual_r1), 1)]
    rm(qual_r1)
    
    # Major ROI enriched columns
    dt_final[, major_roi_method := major_roi_result$method]
    dt_final[, major_roi_matched_seq := major_roi_result$matched_seq]
    dt_final[, `:=`(
      major_roi_n_substitutions = major_roi_result$n_substitutions,
      major_roi_n_deletions = major_roi_result$n_deletions,
      major_roi_n_insertions = major_roi_result$n_insertions
    )]
    
    # Minor ROI enriched columns
    if (use_minor_roi) {
      dt_final[, minor_roi_method := minor_roi_result$method]
      dt_final[, minor_roi_matched_seq := minor_roi_result$matched_seq]
      dt_final[, `:=`(
        minor_roi_n_substitutions = minor_roi_result$n_substitutions,
        minor_roi_n_deletions = minor_roi_result$n_deletions,
        minor_roi_n_insertions = minor_roi_result$n_insertions
      )]
    } else {
      dt_final[, `:=`(
        minor_roi_method = NA_character_,
        minor_roi_matched_seq = NA_character_,
        minor_roi_n_substitutions = NA_integer_,
        minor_roi_n_deletions = NA_integer_,
        minor_roi_n_insertions = NA_integer_
      )]
    }
    
    # ROI Q-scores — vectorized: extract sublists for ROI regions
    # Helper: compute mean qscore for a region defined by start/end vectors
    .roi_qscores_vec <- function(qual_int, qual_lens, starts, ends) {
      n <- length(starts)
      qs <- rep(NA_real_, n)
      valid <- !is.na(starts) & !is.na(ends) & starts >= 1L
      if (!any(valid)) return(qs)
      v_idx <- which(valid)
      capped_ends <- pmin(ends[v_idx], qual_lens[v_idx])
      roi_ranges <- IRanges::IRanges(start = starts[v_idx], end = capped_ends)
      offsets <- c(0L, cumsum(qual_lens)[-length(qual_lens)])
      shifted <- IRanges::shift(roi_ranges, shift = offsets[v_idx])
      flat_qual <- unlist(qual_int, use.names = FALSE)
      roi_quals <- IRanges::extractList(flat_qual, shifted)
      qs[v_idx] <- round(sum(roi_quals) / lengths(roi_quals), 1)
      return(qs)
    }
    
    # Helper: extract single Q-score at a position
    .base_qscores_vec <- function(qual_int, qual_lens, positions) {
      n <- length(positions)
      qs <- rep(NA_real_, n)
      valid <- !is.na(positions) & positions >= 1L & positions <= qual_lens
      if (!any(valid)) return(qs)
      v_idx <- which(valid)
      offsets <- c(0L, cumsum(qual_lens)[-length(qual_lens)])
      flat_pos <- offsets[v_idx] + positions[v_idx]
      flat_qual <- unlist(qual_int, use.names = FALSE)
      qs[v_idx] <- as.numeric(flat_qual[flat_pos])
      return(qs)
    }
    
    dt_final[, major_roi_qscore := .roi_qscores_vec(
      qual_int, qual_lens, major_roi_result$roi_start, major_roi_result$roi_end)]
    
    major_n_read_pos <- major_roi_result$roi_start + major_roi_result$n_offset - 1L
    dt_final[, major_base_qscore := .base_qscores_vec(qual_int, qual_lens, major_n_read_pos)]
    
    if (use_minor_roi) {
      dt_final[, minor_roi_qscore := .roi_qscores_vec(
        qual_int, qual_lens, minor_roi_result$roi_start, minor_roi_result$roi_end)]
      minor_n_read_pos <- minor_roi_result$roi_start + minor_roi_result$n_offset - 1L
      dt_final[, minor_base_qscore := .base_qscores_vec(qual_int, qual_lens, minor_n_read_pos)]
    } else {
      dt_final[, `:=`(minor_roi_qscore = NA_real_, minor_base_qscore = NA_real_)]
    }
    
    rm(qual_bstring, qual_int, qual_lens)
  }
  
  n_major_detected <- sum(!is.na(major_snv_base))
  n_minor_detected <- sum(!is.na(minor_snv_base))
  message(sprintf(
    "  Chunk: %d -> %d len -> %d R1 -> %d chimeric -> %d polyT -> %d SNV (%.0f%%)",
    n_input, n_length_filter, n_valid_read1, n_chimeric,
    n_polyT_pass, n_major_detected, 100 * n_major_detected / n_polyT_pass
  ))
  
  # Column ordering
  core_cols <- c("id", "CBC", "UMI", "strand",
                 "major_base", "major_status", "minor_base", "minor_status")
  if (enriched) {
    enriched_cols <- c("read_length", "read1_pos", "read1_match_width",
                       "read1_n_substitutions", "read1_n_deletions", "read1_n_insertions",
                       "read1_qscore", "mean_qscore",
                       "major_roi_method", "major_roi_matched_seq",
                       "major_roi_n_substitutions", "major_roi_n_deletions", "major_roi_n_insertions",
                       "major_roi_qscore", "major_base_qscore",
                       "minor_roi_method", "minor_roi_matched_seq",
                       "minor_roi_n_substitutions", "minor_roi_n_deletions", "minor_roi_n_insertions",
                       "minor_roi_qscore", "minor_base_qscore")
    data.table::setcolorder(dt_final, c(core_cols, enriched_cols))
  } else {
    data.table::setcolorder(dt_final, core_cols)
  }
  data.table::setnames(dt_final, "id", "read_id")
  return(list(data = dt_final, stats = list(
    n_input = n_input, n_length_filter = n_length_filter, n_valid_read1 = n_valid_read1,
    n_chimeric = n_chimeric, n_polyT_pass = n_polyT_pass,
    n_major_detected = n_major_detected, n_minor_detected = n_minor_detected)))
}


# --- Main function ---

#' Detect SNVs in Nanopore long-read scRNA-seq data
#'
#' Streams through a FASTQ file, detects adapter sequences, extracts cell
#' barcodes and UMIs, and identifies SNVs at user-defined regions of interest
#' using a two-pass strategy (exact match + pairwise alignment).
#'
#' @param session_name Character. Prefix for output files.
#' @param path_input_file Character. Path to FASTQ file (.fastq or .fastq.gz).
#' @param path_output_folder Character. Directory for output files.
#' @param read1_seq Character. Adapter sequence before the barcode.
#' @param read1_mismatch Integer. Allowed mismatches in adapter detection.
#' @param read1_within_range Integer. Search window from read start (default: 70).
#' @param polyN_base Character or NULL. Base of the homopolymer tail required
#'   after the UMI (default: "T", i.e. the poly(T) of 3' GEX libraries). Set to
#'   \code{NULL} to disable the tail filter, e.g. for 5' GEX libraries where the
#'   template-switch oligo rather than a poly(T) follows the UMI.
#' @param polyN_length Integer. Required poly-N length (default: 10). Ignored
#'   when \code{polyN_base = NULL}.
#' @param polyN_mismatch Integer. Allowed mismatches in poly-N (default: 2).
#'   Ignored when \code{polyN_base = NULL}.
#' @param major_roi_seq Character. ROI sequence with 'N' at mutation position.
#' @param major_roi_mismatch Integer. Allowed mismatches for major ROI.
#' @param major_roi_mut Character or NULL. Base representing mutation (e.g. "T").
#' @param major_roi_wt Character or NULL. Base representing wild-type (e.g. "A").
#' @param minor_roi_seq Character or NULL. Secondary ROI for SNP phasing.
#' @param minor_roi_mismatch Integer or NULL. Mismatches for minor ROI.
#' @param minor_roi_mut Character or NULL. Minor ROI mutation base.
#' @param minor_roi_wt Character or NULL. Minor ROI wild-type base.
#' @param len_CBC Integer. Cell barcode length (default: 16).
#' @param len_UMI Integer. UMI length (default: 12).
#' @param n_cores Integer. CPU cores for parallel processing.
#' @param n_reads_per_chunk Integer. Reads per streaming chunk (default: 1e5).
#' @param enriched Logical. If TRUE, output includes additional columns for
#'   publication analyses: read_length, read1_pos, read1_match_width,
#'   read1_n_substitutions, read1_n_deletions, read1_n_insertions,
#'   read1_qscore, mean_qscore, major/minor_roi_method,
#'   major/minor_roi_matched_seq, major/minor_roi_n_substitutions,
#'   major/minor_roi_n_deletions, major/minor_roi_n_insertions,
#'   major/minor_roi_qscore, major/minor_base_qscore (default: FALSE).
#'
#' @return Invisibly returns NULL. Writes a CSV file
#'   \code{{session_name}_snv_table.csv}. Core columns: read_id, CBC, UMI,
#'   strand, major_base, major_status, minor_base, minor_status. When
#'   \code{enriched = TRUE}, appends 22 additional columns with detection
#'   metadata, quality scores, per-error-type counts, and per-region/base
#'   Q-scores.
#'
#' @examples
#' \dontrun{
#' detect_snv(
#'     session_name       = "EXP28_iPSC",
#'     path_input_file    = "experiment.fastq.gz",
#'     path_output_folder = "output/",
#'     read1_seq          = "CTACACGACGCTCTTCCGATCT",
#'     read1_mismatch     = 4L,
#'     major_roi_seq      = "GAGATTTCNCTGTAGCT",
#'     major_roi_mismatch = 3L,
#'     major_roi_mut      = "T",
#'     major_roi_wt       = "A",
#'     n_cores            = 4L
#' )
#' }
#'
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
                       n_reads_per_chunk = 1e5,
                       enriched = FALSE) {
  stopifnot("`session_name` must be character" = is.character(session_name))
  stopifnot("`session_name` must not contain special characters" = grepl("^[a-zA-Z0-9_-]+$", session_name))
  stopifnot("`path_input_file` must be character" = is.character(path_input_file))
  stopifnot("`path_input_file` - file not found" = file.exists(path_input_file))
  stopifnot("`path_input_file` must be .fastq or .fastq.gz" = grepl("\\.fastq(\\.gz)?$", path_input_file))
  stopifnot("`path_output_folder` must be character" = is.character(path_output_folder))
  stopifnot("`path_output_folder` - path not found" = dir.exists(path_output_folder))
  stopifnot("`read1_seq` must contain only A,C,G,T" = grepl("^[GATC]+$", read1_seq))
  stopifnot("`read1_mismatch` must be integer" = is.numeric(read1_mismatch) && read1_mismatch %% 1 == 0)
  stopifnot("`major_roi_seq` must have exactly one N" = grepl("^[GATCRYSWKMBDHV]*N[GATCRYSWKMBDHV]*$", major_roi_seq))
  stopifnot("`major_roi_mismatch` must be integer" = is.numeric(major_roi_mismatch) && major_roi_mismatch %% 1 == 0)
  stopifnot("`n_cores` must be integer" = is.numeric(n_cores) && n_cores %% 1 == 0)
  
  if (is.null(polyN_base)) {
    message("  Poly(N) tail filter: DISABLED (polyN_base = NULL)")
  } else {
    stopifnot("`polyN_base` must be a single base (A, C, G or T) or NULL" = grepl("^[ACGT]$", polyN_base))
    stopifnot("`polyN_length` must be a positive integer" =
                is.numeric(polyN_length) && polyN_length %% 1 == 0 && polyN_length > 0)
    stopifnot("`polyN_mismatch` must be a non-negative integer" =
                is.numeric(polyN_mismatch) && polyN_mismatch %% 1 == 0 && polyN_mismatch >= 0)
  }
  
  use_major_labeling <- !is.null(major_roi_mut) || !is.null(major_roi_wt)
  if (use_major_labeling) {
    stopifnot("Both major_roi_mut and major_roi_wt must be provided together" =
                !is.null(major_roi_mut) && !is.null(major_roi_wt))
    stopifnot("`major_roi_mut` must be single nucleotide" = grepl("^[ACGT]$", major_roi_mut))
    stopifnot("`major_roi_wt` must be single nucleotide" = grepl("^[ACGT]$", major_roi_wt))
    stopifnot("`major_roi_mut` and `major_roi_wt` must differ" = major_roi_mut != major_roi_wt)
    message("  Major ROI labeling: ", major_roi_mut, " = MUT, ", major_roi_wt, " = WT")
  }
  
  use_minor_roi <- !is.null(minor_roi_seq)
  use_minor_labeling <- FALSE
  if (use_minor_roi) {
    stopifnot("`minor_roi_seq` must have exactly one N" = grepl("^[GATCRYSWKMBDHV]*N[GATCRYSWKMBDHV]*$", minor_roi_seq))
    stopifnot("`minor_roi_mismatch` required" = !is.null(minor_roi_mismatch))
    use_minor_labeling <- !is.null(minor_roi_mut) || !is.null(minor_roi_wt)
    if (use_minor_labeling) {
      stopifnot("Both minor_roi_mut and minor_roi_wt must be provided together" =
                  !is.null(minor_roi_mut) && !is.null(minor_roi_wt))
      message("  Minor ROI labeling: ", minor_roi_mut, " = MUT, ", minor_roi_wt, " = WT")
    }
  } else {
    message("  Minor ROI: not provided, skipping")
  }
  
  available_cores <- parallel::detectCores()
  if (n_cores > available_cores - 1) {
    n_cores <- max(1L, available_cores - 1L)
    message("  n_cores adjusted to ", n_cores)
  }
  
  read1_seq <- Biostrings::DNAString(read1_seq)
  major_roi_seq <- Biostrings::DNAString(major_roi_seq)
  if (use_minor_roi) minor_roi_seq <- Biostrings::DNAString(minor_roi_seq)
  
  major_snv_cache <- new.env(hash = TRUE, parent = emptyenv())
  minor_snv_cache <- if (use_minor_roi) new.env(hash = TRUE, parent = emptyenv()) else NULL
  read1_error_cache <- if (enriched) new.env(hash = TRUE, parent = emptyenv()) else NULL
  
  out_file <- file.path(path_output_folder, paste0(session_name, "_snv_table.csv"))
  
  start_time <- Sys.time()
  message(start_time, " - Starting SNV detection")
  message("  Input: ", path_input_file)
  message("  Chunk size: ", format(n_reads_per_chunk, scientific = FALSE))
  if (enriched) message("  Enriched mode: ON (additional columns in output)")
  
  fq_stream <- ShortRead::FastqStreamer(path_input_file, n = n_reads_per_chunk)
  on.exit(close(fq_stream), add = TRUE)
  
  chunk_id <- 1L
  total_reads <- 0L
  total_detected <- 0L
  total_length_filter <- 0L
  total_valid_read1 <- 0L
  total_chimeric <- 0L
  total_polyT_pass <- 0L
  total_major_detected <- 0L
  total_minor_detected <- 0L
  
  repeat {
    fq_chunk <- ShortRead::yield(fq_stream)
    if (length(fq_chunk) == 0L) {
      message("All chunks processed.")
      break
    }
    
    chunk_result <- .process_chunk(
      fq_chunk,
      read1_seq, read1_mismatch, read1_within_range,
      polyN_base, polyN_length, polyN_mismatch,
      major_roi_seq, major_roi_mismatch,
      major_roi_mut, major_roi_wt,
      minor_roi_seq, minor_roi_mismatch,
      minor_roi_mut, minor_roi_wt,
      len_CBC, len_UMI, n_cores, major_snv_cache, minor_snv_cache,
      read1_error_cache,
      use_minor_roi, use_major_labeling, use_minor_labeling,
      enriched
    )
    
    results <- chunk_result$data
    cs <- chunk_result$stats
    total_length_filter <- total_length_filter + cs$n_length_filter
    total_valid_read1 <- total_valid_read1 + cs$n_valid_read1
    total_chimeric <- total_chimeric + cs$n_chimeric
    total_polyT_pass <- total_polyT_pass + cs$n_polyT_pass
    total_major_detected <- total_major_detected + cs$n_major_detected
    total_minor_detected <- total_minor_detected + cs$n_minor_detected
    
    if (nrow(results) > 0) {
      data.table::fwrite(results, out_file, append = (chunk_id > 1L))
      total_detected <- total_detected + nrow(results)
    }
    
    total_reads <- total_reads + length(fq_chunk)
    message(sprintf("%s - Chunk %d: %d reads processed, %d total detected",
                    Sys.time(), chunk_id, length(fq_chunk), total_detected))
    
    chunk_id <- chunk_id + 1L
    gc(verbose = FALSE)
  }
  
  end_time <- Sys.time()
  runtime_secs <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 1)
  
  message("\n--- Complete ---")
  message("Total reads: ", format(total_reads, big.mark = ","))
  message("Length filtered: ", format(total_length_filter, big.mark = ","))
  message("Valid Read1: ", format(total_valid_read1, big.mark = ","))
  message("Chimeric: ", format(total_chimeric, big.mark = ","))
  message("PolyT pass: ", format(total_polyT_pass, big.mark = ","))
  message("Major SNV detected: ", format(total_major_detected, big.mark = ","))
  message("Minor SNV detected: ", format(total_minor_detected, big.mark = ","))
  message("Cache unique (major): ", length(major_snv_cache))
  if (!is.null(minor_snv_cache)) message("Cache unique (minor): ", length(minor_snv_cache))
  if (!is.null(read1_error_cache)) message("Cache unique (read1): ", length(read1_error_cache))
  message("Runtime: ", runtime_secs, "s")
  message("Output: ", out_file)
  
  # Write log CSV
  log_dt <- data.table::data.table(
    session_name = session_name,
    input_file = basename(path_input_file),
    read1_mismatch = as.integer(read1_mismatch),
    major_roi_mismatch = as.integer(major_roi_mismatch),
    minor_roi_mismatch = if (use_minor_roi) as.integer(minor_roi_mismatch) else NA_integer_,
    enriched = enriched,
    n_cores = as.integer(n_cores),
    total_reads = total_reads,
    total_length_filtered = total_length_filter,
    total_valid_read1 = total_valid_read1,
    total_chimeric = total_chimeric,
    total_polyT_pass = total_polyT_pass,
    total_major_detected = total_major_detected,
    total_minor_detected = total_minor_detected,
    cache_unique_major = length(major_snv_cache),
    cache_unique_minor = if (!is.null(minor_snv_cache)) length(minor_snv_cache) else NA_integer_,
    cache_unique_read1 = if (!is.null(read1_error_cache)) length(read1_error_cache) else NA_integer_,
    runtime_seconds = runtime_secs
  )
  log_file <- file.path(path_output_folder, paste0(session_name, "_log.csv"))
  data.table::fwrite(log_dt, log_file)
  message("Log: ", log_file)
  
  invisible(NULL)
}
