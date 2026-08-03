# =============================================================================
# flag_snv.R - SNV mutation flagging with total-read thresholds (MERLIN)
# =============================================================================
#
# Classifies each cell barcode as MUT / WT / MUT_WT from the per-cell UMI summary
# produced by summarize_snv(). UMIs are first filtered by their total read count,
# then flagged MUT/WT by a strict consensus, and finally combined into a cell
# call. When a minor SNP ROI is available, an optional phasing layer separates
# definitive from ambiguous wild-type cells and validates MUT calls against the
# expected minor allele.
# =============================================================================


# --- Internal helpers --------------------------------------------------------

#' Flag one cell's mutation status at a given read threshold
#'
#' Runs the four-level calling strategy for a single cell: (1) filter UMIs by
#' total reads and flag each as MUT/WT by consensus, (2) optionally phase each
#' UMI against the minor SNP, (3) derive the major cell call, and (4) derive the
#' phased cell call.
#'
#' @param mut_calls_str,wt_calls_str Comma-separated per-UMI MUT/WT read counts.
#' @param snp_mut_str,snp_wt_str Comma-separated per-UMI minor-ROI counts
#'   (only used when \code{use_phasing = TRUE}).
#' @param threshold Minimum total reads (mut + wt) for a UMI to be considered.
#' @param consensus_threshold Fraction of reads a UMI must lean to be called.
#' @param use_phasing Whether to apply the minor-SNP phasing logic.
#' @param minor_percentage Minimum minor reads relative to major reads.
#' @param minor_on_major Expected minor allele on the major allele ("MUT"/"WT").
#' @return List with elements \code{major} and \code{phased}.
#' @noRd
.flag_cell_at_threshold <- function(mut_calls_str, wt_calls_str,
                                    snp_mut_str = NULL, snp_wt_str = NULL,
                                    threshold,
                                    consensus_threshold = 0.8,
                                    use_phasing = FALSE,
                                    minor_percentage = 0.5,
                                    minor_on_major = "MUT") {
    # Parse the per-UMI read counts from their comma-separated strings.
    mut_vec <- as.numeric(str_split(mut_calls_str, ",", simplify = TRUE))
    wt_vec <- as.numeric(str_split(wt_calls_str, ",", simplify = TRUE))

    # No usable data for this cell.
    if (length(mut_vec) == 0 || all(is.na(mut_vec))) {
        return(list(major = "no_data", phased = "no_data"))
    }

    mut_vec[is.na(mut_vec)] <- 0
    wt_vec[is.na(wt_vec)] <- 0

    total_per_umi <- mut_vec + wt_vec

    # -------------------------------------------------------------------------
    # LEVEL 1: filter UMIs by total reads, then flag each MUT/WT by consensus
    # -------------------------------------------------------------------------

    # Every UMI (a unique cbc+umi) must clear the total-read threshold.
    valid_umi_idx <- total_per_umi >= threshold

    if (!any(valid_umi_idx)) {
        return(list(major = "below_threshold", phased = "below_threshold"))
    }

    mut_vec <- mut_vec[valid_umi_idx]
    wt_vec <- wt_vec[valid_umi_idx]
    total_per_umi <- total_per_umi[valid_umi_idx]

    # A UMI is called only if its reads lean >= consensus_threshold one way.
    pct_mut <- mut_vec / total_per_umi
    pct_wt <- wt_vec / total_per_umi

    umi_calls <- rep("ambiguous", length(total_per_umi))

    umi_calls[pct_mut >= consensus_threshold] <- "MUT"
    umi_calls[pct_wt >= consensus_threshold] <- "WT"

    # Only decisively-called UMIs contribute to the cell status.
    valid_calls <- umi_calls[umi_calls %in% c("MUT", "WT")]

    # Positions of those decisive calls within the threshold-passing UMIs.
    valid_call_indices <- which(umi_calls %in% c("MUT", "WT"))

    if (length(valid_calls) == 0) {
        return(list(major = "below_threshold", phased = "below_threshold"))
    }

    # -------------------------------------------------------------------------
    # LEVEL 2: per-UMI phasing against the minor SNP
    # -------------------------------------------------------------------------

    phased_calls <- rep("ambiguous", length(valid_calls))

    if (use_phasing && !is.null(snp_mut_str) && !is.null(snp_wt_str)) {
        snp_mut <- as.numeric(str_split(snp_mut_str, ",", simplify = TRUE))
        snp_wt <- as.numeric(str_split(snp_wt_str, ",", simplify = TRUE))

        snp_mut[is.na(snp_mut)] <- 0
        snp_wt[is.na(snp_wt)] <- 0

        # Align the minor counts to the same threshold-passing UMIs...
        if (length(snp_mut) >= length(valid_umi_idx)) {
            snp_mut <- snp_mut[valid_umi_idx]
            snp_wt <- snp_wt[valid_umi_idx]
        }

        # ...then to only those UMIs that received a decisive MUT/WT call.
        snp_mut <- snp_mut[valid_call_indices]
        snp_wt <- snp_wt[valid_call_indices]
        total_minor <- snp_mut + snp_wt

        # Matching major totals, for the minor-vs-major percentage check.
        major_totals <- total_per_umi[valid_call_indices]

        # Require enough minor reads relative to the major signal to trust it.
        has_enough_minor <- total_minor >= (major_totals * minor_percentage)

        # Consensus direction of the minor SNP within each UMI.
        pct_minor_mut <- snp_mut / total_minor
        pct_minor_wt <- snp_wt / total_minor

        # Combine major call x minor call into a per-UMI phased label.
        for (i in seq_along(valid_calls)) {
            major_status <- valid_calls[i]

            # "none" = insufficient or non-consensus minor signal.
            minor_status <- "none"
            if (has_enough_minor[i]) {
                if (!is.na(pct_minor_mut[i]) && pct_minor_mut[i] >= consensus_threshold) {
                    minor_status <- "MUT"
                } else if (!is.na(pct_minor_wt[i]) && pct_minor_wt[i] >= consensus_threshold) {
                    minor_status <- "WT"
                }
            }

            # The minor allele expected to sit on the major (mutant) allele.
            expected_minor <- minor_on_major
            other_minor <- ifelse(expected_minor == "MUT", "WT", "MUT")

            if (major_status == "MUT") {
                if (minor_status == "none") {
                    phased_calls[i] <- "call_MUT_non_mnsnp" # MUT, minor SNP not seen
                } else if (minor_status == expected_minor) {
                    phased_calls[i] <- "call_MUT_mnsnp" # MUT confirmed by minor SNP
                } else if (minor_status == other_minor) {
                    phased_calls[i] <- "call_MUT_ERROR" # MUT contradicted by minor SNP
                }
            } else if (major_status == "WT") {
                if (minor_status == "none") {
                    phased_calls[i] <- "amb_WT" # WT, cannot be phased -> ambiguous
                } else if (minor_status == expected_minor) {
                    phased_calls[i] <- "def_WT" # WT confirmed by minor SNP
                } else if (minor_status == other_minor) {
                    phased_calls[i] <- "amb_WT"
                }
            }
        }
    } else {
        # Without phasing data, MUT UMIs are unverified and WT UMIs are ambiguous.
        phased_calls <- ifelse(valid_calls == "MUT", "call_MUT_non_mnsnp", "amb_WT")
    }

    # -------------------------------------------------------------------------
    # LEVEL 3: major cell call from the decisive UMI calls
    # -------------------------------------------------------------------------

    has_mut <- any(valid_calls == "MUT")
    has_wt <- any(valid_calls == "WT")

    major_final <- "below_threshold"
    if (has_mut && has_wt) {
        major_final <- "MUT_WT" # heterozygous: both alleles present
    } else if (has_mut) {
        major_final <- "MUT"
    } else if (has_wt) {
        major_final <- "WT"
    }

    # -------------------------------------------------------------------------
    # LEVEL 4: phased cell call, aggregating the per-UMI phased labels
    # -------------------------------------------------------------------------

    phased_final <- "amb_WT" # default

    count_mut_error <- sum(phased_calls == "call_MUT_ERROR")
    count_mut_mnsnp <- sum(phased_calls == "call_MUT_mnsnp")
    count_mut_non_mnsnp <- sum(phased_calls == "call_MUT_non_mnsnp")

    # Older WT-branch flags, kept for compatibility with mixed inputs.
    has_error_wt <- any(phased_calls == "ERROR")
    has_def_wt <- any(phased_calls == "def_WT")

    # Any (verified or unverified) MUT evidence at the UMI level.
    has_any_mut_call <- (count_mut_mnsnp > 0 || count_mut_non_mnsnp > 0)

    if (major_final == "MUT") {
        # Pure-MUT cell: resolve which kind of MUT evidence it carries.
        if (count_mut_error > 0) {
            phased_final <- "MUT_ERROR"
        } else if (count_mut_mnsnp > 0 && count_mut_non_mnsnp > 0) {
            phased_final <- "MUT_both"
        } else if (count_mut_mnsnp > 0) {
            phased_final <- "MUT_mnsnp"
        } else if (count_mut_non_mnsnp > 0) {
            phased_final <- "MUT_non_mnsnp"
        } else {
            phased_final <- "MUT_non_mnsnp" # unreachable when major == MUT
        }
    } else {
        # WT or heterozygous (MUT_WT) cell.
        is_error <- has_error_wt || (count_mut_error > 0)

        if (is_error) {
            phased_final <- "ERROR1"
        } else if (has_any_mut_call) {
            # Surface MUT evidence even in heterozygous / ambiguous-WT cells.
            if (count_mut_mnsnp > 0 && count_mut_non_mnsnp > 0) {
                phased_final <- "MUT_both"
            } else if (count_mut_mnsnp > 0) {
                phased_final <- "MUT_mnsnp"
            } else {
                phased_final <- "MUT_non_mnsnp"
            }
        } else if (has_def_wt) {
            phased_final <- "def_WT"
        } else {
            phased_final <- "amb_WT"
        }
    }

    if (length(phased_calls) == 0) phased_final <- "below_threshold"

    return(list(major = major_final, phased = phased_final))
}


# --- Main function -----------------------------------------------------------

#' Flag SNV mutation status per cell at a chosen threshold
#'
#' Classifies each cell barcode from the per-cell UMI summary produced by
#' \code{\link{summarize_snv}}. UMIs are filtered by total read count, flagged
#' MUT/WT by strict consensus, and aggregated into a major call
#' (\code{MUT}/\code{WT}/\code{MUT_WT}). With \code{use_phasing = TRUE}, a minor
#' SNP ROI is used to split wild-type cells into definitive vs. ambiguous and to
#' validate mutant calls.
#'
#' @param path_snv_summary Path to the summary CSV from \code{\link{summarize_snv}},
#'   or the \code{data.table} returned as \code{summarize_snv()$summary}.
#' @param path_barcodes Optional full barcode whitelist used to backfill missing
#'   cells so every expected barcode appears in the output.
#' @param path_barcode_filter Optional whitelist; when supplied, ONLY these
#'   barcodes are analysed.
#' @param path_output_folder Optional output directory; if \code{NULL}, nothing
#'   is written to disk.
#' @param session_name Prefix for output files (default: \code{"snv_flags"}).
#' @param threshold Minimum total reads per UMI to consider (default: 10).
#' @param consensus_threshold Fraction of reads a UMI must lean to be called
#'   (default: 0.8).
#' @param use_phasing Use minor-SNP data for definitive/ambiguous calls
#'   (default: FALSE).
#' @param minor_percentage Minimum fraction of minor reads relative to major
#'   reads (default: 0.2).
#' @param minor_on_major Expected minor allele on the major allele, \code{"MUT"}
#'   or \code{"WT"} (default: \code{"MUT"}).
#'
#' @return A \code{data.table} with one row per cell barcode and columns
#'   \code{cbc}, \code{barcode_10x}, \code{has_np_data}, \code{mutation_call},
#'   \code{major_call_simple}, \code{phased_call} and \code{phased_call_simple}.
#'
#' @seealso \code{\link{detect_snv}}, \code{\link{summarize_snv}}
#' @export
flag_snv <- function(path_snv_summary,
                     path_barcodes = NULL,
                     path_barcode_filter = NULL,
                     path_output_folder = NULL,
                     session_name = "snv_flags",
                     threshold = 10,
                     consensus_threshold = 0.8,
                     use_phasing = FALSE,
                     minor_percentage = 0.2,
                     minor_on_major = "MUT") {
    # --- Load the per-cell summary (file path or in-memory data.table) ---
    message(Sys.time(), " - Loading SNV summary for flagging")

    if (is.character(path_snv_summary)) {
        stopifnot("Summary file not found" = file.exists(path_snv_summary))
        dt_summary <- fread(path_snv_summary, sep = ";")
    } else if (is.data.table(path_snv_summary)) {
        dt_summary <- copy(path_snv_summary)
    } else {
        stop("path_snv_summary must be file path or data.table")
    }

    req_cols <- c("cbc", "n_umi_mut", "n_umi_wt")
    missing <- setdiff(req_cols, colnames(dt_summary))
    if (length(missing) > 0) {
        stop("Missing required columns: ", paste(missing, collapse = ", "))
    }

    message("  Loaded ", nrow(dt_summary), " cell barcodes")

    # Phasing needs the minor-SNP columns; disable it (with a warning) if absent.
    has_snp_data <- all(c("n_snp_mut", "n_snp_wt") %in% colnames(dt_summary))
    if (use_phasing && !has_snp_data) {
        warning("use_phasing=TRUE but no SNP columns found, proceeding without phasing")
        use_phasing <- FALSE
    }

    # --- Resolve the set of barcodes to report (backfill and/or filter) ---

    # 1. Start from the barcodes present in the summary.
    all_barcodes <- dt_summary$cbc

    # 2. Optionally add a backfill list so ungenotyped cells still appear.
    if (!is.null(path_barcodes) && file.exists(path_barcodes)) {
        message("  Loading full barcode list (backfill): ", path_barcodes)
        if (grepl("\\.csv$", path_barcodes, ignore.case = TRUE)) {
            bc_raw <- fread(path_barcodes, header = FALSE)[[1]]
        } else {
            bc_raw <- readLines(path_barcodes)
        }
        all_barcodes <- unique(c(substr(bc_raw, 1, 16), all_barcodes))
    }

    # 3. Optionally restrict to a filter list (intersection).
    if (!is.null(path_barcode_filter) && file.exists(path_barcode_filter)) {
        message("  Loading barcode filter list: ", path_barcode_filter)
        if (grepl("\\.csv$", path_barcode_filter, ignore.case = TRUE)) {
            bc_filt_raw <- fread(path_barcode_filter, header = FALSE)[[1]]
        } else {
            bc_filt_raw <- readLines(path_barcode_filter)
        }
        filter_set <- unique(substr(bc_filt_raw, 1, 16))

        original_count <- length(all_barcodes)
        all_barcodes <- intersect(all_barcodes, filter_set)

        message(sprintf("  Filtered barcodes: %d -> %d", original_count, length(all_barcodes)))

        # Drop unwanted cells from the summary too, so they are not processed.
        dt_summary <- dt_summary[cbc %in% all_barcodes]
    }

    message("  Total unique barcodes to process: ", length(all_barcodes))

    # --- Flag mutations for every barcode ---
    message(Sys.time(), " - Flagging mutations")
    message("  Total read threshold: ", threshold)
    message("  Consensus threshold: ", consensus_threshold * 100, "%")
    if (use_phasing) {
        message("  Phasing: ENABLED")
        message("  Minor read percentage: ", minor_percentage * 100, "%")
        message("  Minor on major allele: ", minor_on_major)
    }

    # One row per barcode, joined to its summary (missing cells -> NA columns).
    dt_calls <- data.table(cbc = all_barcodes)

    setkey(dt_calls, cbc)
    setkey(dt_summary, cbc)
    dt_calls <- merge(dt_calls, dt_summary, all.x = TRUE)

    # Cells with no Nanopore data have no MUT counts.
    dt_calls[, has_np_data := !is.na(n_umi_mut)]

    # Compute the major and phased calls, row by row.
    dt_calls[, c("mutation_call", "phased_call") := {
        if (!has_np_data) {
            list("no_np_data", "no_np_data")
        } else {
            res_list <- mapply(
                .flag_cell_at_threshold,
                n_umi_mut, n_umi_wt,
                if (use_phasing) n_snp_mut else NA,
                if (use_phasing) n_snp_wt else NA,
                MoreArgs = list(
                    threshold = threshold,
                    consensus_threshold = consensus_threshold,
                    use_phasing = use_phasing,
                    minor_percentage = minor_percentage,
                    minor_on_major = minor_on_major
                ),
                SIMPLIFY = FALSE
            )
            # Split the list-of-lists into the two output columns.
            list(
                sapply(res_list, `[[`, "major"),
                sapply(res_list, `[[`, "phased")
            )
        }
    }, by = seq_len(nrow(dt_calls))]

    # Log the major call distribution.
    message("  Major Call Distribution:")
    call_counts <- table(dt_calls$mutation_call)
    message(sprintf(
        "    MUT: %d | WT: %d | MUT_WT: %d | below: %d",
        sum(call_counts["MUT"], na.rm = TRUE),
        sum(call_counts["WT"], na.rm = TRUE),
        sum(call_counts["MUT_WT"], na.rm = TRUE),
        sum(call_counts["below_threshold"], na.rm = TRUE)
    ))

    if (use_phasing) {
        message("  Phased Call Distribution:")
        phased_counts <- table(dt_calls$phased_call)
        print(phased_counts)
    }

    # --- Simplified call columns for convenience ---

    # Major: collapse MUT and heterozygous MUT_WT into a single "MUT".
    dt_calls[, major_call_simple := NA_character_]
    dt_calls[mutation_call == "WT", major_call_simple := "WT"]
    dt_calls[mutation_call %in% c("MUT", "MUT_WT"), major_call_simple := "MUT"]

    # Phased: collapse the MUT_* labels into MUT and map WT confidence levels.
    dt_calls[, phased_call_simple := NA_character_]
    dt_calls[phased_call %in% c("MUT_mnsnp", "MUT_non_mnsnp", "MUT_both"), phased_call_simple := "MUT"]
    dt_calls[phased_call == "amb_WT", phased_call_simple := "ambWT"]
    dt_calls[phased_call == "def_WT", phased_call_simple := "defWT"]

    # Select and order the output columns.
    dt_output <- dt_calls[, .(cbc, has_np_data, mutation_call, major_call_simple, phased_call, phased_call_simple)]

    # Re-attach the 10x-style barcode suffix for downstream joins.
    dt_output[, barcode_10x := paste0(cbc, "-1")]
    setcolorder(dt_output, c("cbc", "barcode_10x", "has_np_data", "mutation_call", "major_call_simple", "phased_call", "phased_call_simple"))

    # --- Write output ---
    if (!is.null(path_output_folder)) {
        stopifnot("Output folder not found" = dir.exists(path_output_folder))

        message(Sys.time(), " - Saving outputs")

        out_calls <- file.path(path_output_folder, paste0(session_name, "_flags.csv"))
        fwrite(dt_output, out_calls)
        message("  Saved: ", out_calls)
    }

    return(dt_output)
}
