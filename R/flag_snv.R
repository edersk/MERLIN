# =============================================================================
# flag_snv.R - SNV Muation Flagging with Total Read Thresholds
# =============================================================================
#
# Description:
#   Flags mutation status (MUT/WT) for each cell barcode based on UMI summary
#   from summarize_snv(). Uses total read count filtering per UMI and
#   strict consensus thresholds.
#
# Author: Sebastian Eder <sebastian.eder@example.com>
# Created: 2026-01-30
# Version: 1.1.0
#
# Pipeline:
#   detect_snv() -> summarize_snv() -> flag_snv()
#
# Features:
#   - Filter UMIs by TOTAL read count (mut + wt) >= threshold
#   - Flag UMIs as MUT/WT based on strict consensus (lean 80%)
#   - Cell-level calling: MUT, WT, or MUT_WT (heterozygous)
#   - Supports detailed SNP phasing integration:
#       - Minor percentage threshold
#       - Minor on major allele validation
#       - Hierarchical cell calling (ERROR1, ERROR2, MUT, def_WT, amb_WT)
#
# License: MIT License
#   Copyright (c) 2026 Sebastian Eder#
# Dependencies:
#   - data.table (>= 1.14.0)
#   - stringr (>= 1.4.0)
#
# =============================================================================


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Flags mutation status for a single cell at a given threshold
.flag_cell_at_threshold <- function(mut_calls_str, wt_calls_str,
                                    snp_mut_str = NULL, snp_wt_str = NULL,
                                    threshold,
                                    consensus_threshold = 0.8,
                                    use_phasing = FALSE,
                                    minor_percentage = 0.5,
                                    minor_on_major = "MUT") {
    # Parse comma-separated read counts per UMI
    mut_vec <- as.numeric(str_split(mut_calls_str, ",", simplify = TRUE))
    wt_vec <- as.numeric(str_split(wt_calls_str, ",", simplify = TRUE))

    # Handle empty or NA
    if (length(mut_vec) == 0 || all(is.na(mut_vec))) {
        return(list(major = "no_data", phased = "no_data"))
    }

    mut_vec[is.na(mut_vec)] <- 0
    wt_vec[is.na(wt_vec)] <- 0

    total_per_umi <- mut_vec + wt_vec

    # ===========================================================================
    # LEVEL 1: Filter & Call UMIs
    # ===========================================================================

    # Filter UMIs by TOTAL reads threshold
    # "Every major SNV (unique cbc+umi combination) should be filtered for the threshold"
    valid_umi_idx <- total_per_umi >= threshold

    if (!any(valid_umi_idx)) {
        return(list(major = "below_threshold", phased = "below_threshold"))
    }

    mut_vec <- mut_vec[valid_umi_idx]
    wt_vec <- wt_vec[valid_umi_idx]
    total_per_umi <- total_per_umi[valid_umi_idx]

    # Flag UMIs based on direction lean
    # "counts lean 80% (input peramter) in one direction"
    pct_mut <- mut_vec / total_per_umi
    pct_wt <- wt_vec / total_per_umi

    # Initialize calls as ambiguous/uncalled
    umi_calls <- rep("ambiguous", length(total_per_umi))

    # Assign MUT/WT if they pass consensus threshold
    umi_calls[pct_mut >= consensus_threshold] <- "MUT"
    umi_calls[pct_wt >= consensus_threshold] <- "WT"

    # We only care about definite MUT/WT calls for the cell status
    valid_calls <- umi_calls[umi_calls %in% c("MUT", "WT")]

    # Indices of valid calls within the originally valid UMIs
    valid_call_indices <- which(umi_calls %in% c("MUT", "WT"))

    if (length(valid_calls) == 0) {
        return(list(major = "below_threshold", phased = "below_threshold"))
    }

    # ===========================================================================
    # LEVEL 2: Phasing Logic (Per UMI)
    # ===========================================================================

    phased_calls <- rep("ambiguous", length(valid_calls))

    if (use_phasing && !is.null(snp_mut_str) && !is.null(snp_wt_str)) {
        snp_mut <- as.numeric(str_split(snp_mut_str, ",", simplify = TRUE))
        snp_wt <- as.numeric(str_split(snp_wt_str, ",", simplify = TRUE))

        snp_mut[is.na(snp_mut)] <- 0
        snp_wt[is.na(snp_wt)] <- 0

        # Filter to valid UMIs first
        if (length(snp_mut) >= length(valid_umi_idx)) {
            snp_mut <- snp_mut[valid_umi_idx]
            snp_wt <- snp_wt[valid_umi_idx]
        }

        # Now filter to only UMIs that got a valid MUT/WT call
        snp_mut <- snp_mut[valid_call_indices]
        snp_wt <- snp_wt[valid_call_indices]
        total_minor <- snp_mut + snp_wt

        # Corresponding major totals for percentage check
        major_totals <- total_per_umi[valid_call_indices]

        # Check percentage threshold
        has_enough_minor <- total_minor >= (major_totals * minor_percentage)

        # Determine minor call (MUT/WT) based on consensus
        pct_minor_mut <- snp_mut / total_minor
        pct_minor_wt <- snp_wt / total_minor

        # Loop through valid calls to assign phased status
        for (i in seq_along(valid_calls)) {
            major_status <- valid_calls[i]

            # Use 'none' if below percentage threshold
            minor_status <- "none"
            if (has_enough_minor[i]) {
                if (!is.na(pct_minor_mut[i]) && pct_minor_mut[i] >= consensus_threshold) {
                    minor_status <- "MUT"
                } else if (!is.na(pct_minor_wt[i]) && pct_minor_wt[i] >= consensus_threshold) {
                    minor_status <- "WT"
                }
            }

            # Apply user logic table
            # "minor_snv_status_on_major_snv_allele" (e.g., MUT)
            expected_minor <- minor_on_major
            other_minor <- ifelse(expected_minor == "MUT", "WT", "MUT")

            if (major_status == "MUT") {
                if (minor_status == "none") {
                    phased_calls[i] <- "call_MUT_non_mnsnp"
                } else if (minor_status == expected_minor) {
                    phased_calls[i] <- "call_MUT_mnsnp"
                } else if (minor_status == other_minor) {
                    phased_calls[i] <- "call_MUT_ERROR"
                }
            } else if (major_status == "WT") {
                if (minor_status == "none") {
                    phased_calls[i] <- "amb_WT"
                } else if (minor_status == expected_minor) {
                    phased_calls[i] <- "def_WT"
                } else if (minor_status == other_minor) {
                    phased_calls[i] <- "amb_WT"
                }
            }
        }
    } else {
        # Fallback without phasing data
        # If no phasing, we consider it "none" (unverified)
        phased_calls <- ifelse(valid_calls == "MUT", "call_MUT_non_mnsnp", "amb_WT")
    }

    # ===========================================================================
    # LEVEL 3: Cell-level calling (Major)
    # ===========================================================================

    has_mut <- any(valid_calls == "MUT")
    has_wt <- any(valid_calls == "WT")

    major_final <- "below_threshold"
    if (has_mut && has_wt) {
        major_final <- "MUT_WT"
    } else if (has_mut) {
        major_final <- "MUT"
    } else if (has_wt) {
        major_final <- "WT"
    }

    # ===========================================================================
    # LEVEL 4: Cell-level calling (Phased)
    # ===========================================================================

    phased_final <- "amb_WT" # Default

    # Parse UMI-level flags
    count_mut_error <- sum(phased_calls == "call_MUT_ERROR")
    count_mut_mnsnp <- sum(phased_calls == "call_MUT_mnsnp")
    count_mut_non_mnsnp <- sum(phased_calls == "call_MUT_non_mnsnp")

    # Old flags for backward compatibility/other branches
    has_error_wt <- any(phased_calls == "ERROR") # From WT branch if set
    has_def_wt <- any(phased_calls == "def_WT")

    # "MUT" in old logic meant verified/unverified MUT.
    # Here mapped to any of the call_MUT_* (except maybe error?)
    has_any_mut_call <- (count_mut_mnsnp > 0 || count_mut_non_mnsnp > 0)

    if (major_final == "MUT") {
        # Specific user logic for Major=MUT
        if (count_mut_error > 0) {
            phased_final <- "MUT_ERROR"
        } else if (count_mut_mnsnp > 0 && count_mut_non_mnsnp > 0) {
            phased_final <- "MUT_both"
        } else if (count_mut_mnsnp > 0) {
            phased_final <- "MUT_mnsnp"
        } else if (count_mut_non_mnsnp > 0) {
            phased_final <- "MUT_non_mnsnp"
        } else {
            phased_final <- "MUT_non_mnsnp" # Should not be reached if major=MUT
        }
    } else {
        # Fallback logic for WT or MUT_WT (Heterozygous)
        # Reconstruct ERROR conditions
        is_error <- has_error_wt || (count_mut_error > 0)

        if (is_error) {
            phased_final <- "ERROR1"
        } else if (has_any_mut_call) {
            # Resolve specific MUT status even for non-primary MUT cells (e.g. heterozygous with ambiguous WT)
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


# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Flag SNV Mutation Status from UMI Summary
#'
#' Flags mutation status (MUT/WT) for each cell barcode based on UMI summary
#' from \code{\link{summarize_snv}}. Uses total read count filtering per UMI
#' and strict consensus thresholds. Optionally supports SNP phasing for
#' definitive vs. ambiguous wild-type classification.
#'
#' The function applies a multi-level calling strategy:
#' \enumerate{
#'   \item Filter UMIs by total read count (mut + wt) >= threshold
#'   \item Flag UMIs as MUT/WT based on consensus (default 80\%)
#'   \item Assign cell-level major call: MUT, WT, or MUT_WT (heterozygous)
#'   \item Optionally assign phased call using minor SNP data
#' }
#'
#' @param path_snv_summary Path to SNV summary CSV from \code{summarize_snv()} OR
#'   a \code{data.table} returned by \code{summarize_snv()$summary}.
#' @param path_barcodes Path to full barcode whitelist (optional). Used for
#'   backfilling missing cells so that all expected barcodes appear in output.
#' @param path_barcode_filter Path to barcode whitelist for filtering (optional).
#'   When provided, ONLY these barcodes will be analyzed.
#' @param path_output_folder Path for output files. If \code{NULL}, no files
#'   are written.
#' @param session_name Prefix for output files (default: \code{"snv_flags"}).
#' @param threshold Minimum TOTAL reads per UMI to consider (default: 10).
#' @param consensus_threshold Fraction of reads for UMI consensus (default: 0.8).
#'   UMIs must lean >= this fraction to one side to be called MUT or WT.
#' @param use_phasing Logical; use SNP data for definitive/ambiguous calls
#'   (default: \code{FALSE}).
#' @param minor_percentage Minimum fraction of minor reads relative to major reads
#'   (default: 0.2).
#' @param minor_on_major Expected minor allele status on major allele
#'   (\code{"MUT"} or \code{"WT"}, default: \code{"MUT"}).
#'
#' @return A \code{data.table} with columns:
#'   \describe{
#'     \item{cbc}{Cell barcode (16-mer)}
#'     \item{barcode_10x}{Cell barcode with \code{-1} suffix for 10x compatibility}
#'     \item{has_np_data}{Logical; whether the cell had Nanopore data}
#'     \item{mutation_call}{Major mutation call: MUT, WT, MUT_WT, below_threshold, or no_np_data}
#'     \item{major_call_simple}{Simplified major call: MUT, WT, or NA}
#'     \item{phased_call}{Detailed phased call (e.g., MUT_mnsnp, def_WT, amb_WT)}
#'     \item{phased_call_simple}{Simplified phased call: MUT, defWT, ambWT, or NA}
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage with file path
#' flags <- flag_snv(
#'   path_snv_summary = "snv_summary.csv",
#'   threshold = 10,
#'   consensus_threshold = 0.8
#' )
#'
#' # With phasing enabled
#' flags <- flag_snv(
#'   path_snv_summary = "snv_summary.csv",
#'   path_barcodes = "barcodes.csv",
#'   path_output_folder = "output/",
#'   session_name = "experiment1",
#'   threshold = 10,
#'   use_phasing = TRUE,
#'   minor_percentage = 0.2,
#'   minor_on_major = "MUT"
#' )
#'
#' # Using data.table input from summarize_snv()
#' summary_result <- summarize_snv(...)
#' flags <- flag_snv(
#'   path_snv_summary = summary_result$summary,
#'   threshold = 5
#' )
#' }
#'
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
    # =========================================================================
    # LOAD DATA
    # =========================================================================

    message(Sys.time(), " - Loading SNV summary for flagging")

    # Accept either file path or data.table
    if (is.character(path_snv_summary)) {
        stopifnot("Summary file not found" = file.exists(path_snv_summary))
        dt_summary <- fread(path_snv_summary, sep = ";")
    } else if (is.data.table(path_snv_summary)) {
        dt_summary <- copy(path_snv_summary)
    } else {
        stop("path_snv_summary must be file path or data.table")
    }

    # Validate columns
    req_cols <- c("cbc", "n_umi_mut", "n_umi_wt")
    missing <- setdiff(req_cols, colnames(dt_summary))
    if (length(missing) > 0) {
        stop("Missing required columns: ", paste(missing, collapse = ", "))
    }

    message("  Loaded ", nrow(dt_summary), " cell barcodes")

    # Check for phasing data
    has_snp_data <- all(c("n_snp_mut", "n_snp_wt") %in% colnames(dt_summary))
    if (use_phasing && !has_snp_data) {
        warning("use_phasing=TRUE but no SNP columns found, proceeding without phasing")
        use_phasing <- FALSE
    }

    # =========================================================================
    # LOAD FULL BARCODE LIST (optional)
    # =========================================================================

    # =========================================================================
    # LOAD BARCODE LISTS (Backfill & Filter)
    # =========================================================================

    # 1. Start with barcodes from summary
    all_barcodes <- dt_summary$cbc

    # 2. Add backfill list (optional)
    if (!is.null(path_barcodes) && file.exists(path_barcodes)) {
        message("  Loading full barcode list (backfill): ", path_barcodes)
        if (grepl("\\.csv$", path_barcodes, ignore.case = TRUE)) {
            bc_raw <- fread(path_barcodes, header = FALSE)[[1]]
        } else {
            bc_raw <- readLines(path_barcodes)
        }
        # Union
        all_barcodes <- unique(c(substr(bc_raw, 1, 16), all_barcodes))
    }

    # 3. Apply Filter (optional)
    if (!is.null(path_barcode_filter) && file.exists(path_barcode_filter)) {
        message("  Loading barcode filter list: ", path_barcode_filter)
        if (grepl("\\.csv$", path_barcode_filter, ignore.case = TRUE)) {
            bc_filt_raw <- fread(path_barcode_filter, header = FALSE)[[1]]
        } else {
            bc_filt_raw <- readLines(path_barcode_filter)
        }
        filter_set <- unique(substr(bc_filt_raw, 1, 16))

        # Intersect
        original_count <- length(all_barcodes)
        all_barcodes <- intersect(all_barcodes, filter_set)

        message(sprintf("  Filtered barcodes: %d -> %d", original_count, length(all_barcodes)))

        # Also filter the summary data to avoid processing unwanted cells
        dt_summary <- dt_summary[cbc %in% all_barcodes]
    }

    message("  Total unique barcodes to process: ", length(all_barcodes))

    # =========================================================================
    # FLAG MUTATIONS
    # =========================================================================

    message(Sys.time(), " - Flagging mutations")
    message("  Total read threshold: ", threshold)
    message("  Consensus threshold: ", consensus_threshold * 100, "%")
    if (use_phasing) {
        message("  Phasing: ENABLED")
        message("  Minor read percentage: ", minor_percentage * 100, "%")
        message("  Minor on major allele: ", minor_on_major)
    }

    # Initialize result table
    dt_calls <- data.table(cbc = all_barcodes)

    # Merge with summary data
    setkey(dt_calls, cbc)
    setkey(dt_summary, cbc)
    dt_calls <- merge(dt_calls, dt_summary, all.x = TRUE)

    # Mark cells with no nanopore data
    dt_calls[, has_np_data := !is.na(n_umi_mut)]

    # Compute flags

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
            # Unpack list results
            list(
                sapply(res_list, `[[`, "major"),
                sapply(res_list, `[[`, "phased")
            )
        }
    }, by = seq_len(nrow(dt_calls))]

    # Log call distribution
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

    # =========================================================================
    # CREATE OUTPUT
    # =========================================================================

    # Generate simplified columns
    dt_calls[, major_call_simple := NA_character_]
    dt_calls[mutation_call == "WT", major_call_simple := "WT"]
    dt_calls[mutation_call %in% c("MUT", "MUT_WT"), major_call_simple := "MUT"]

    dt_calls[, phased_call_simple := NA_character_]
    dt_calls[phased_call %in% c("MUT_mnsnp", "MUT_non_mnsnp", "MUT_both"), phased_call_simple := "MUT"]
    dt_calls[phased_call == "amb_WT", phased_call_simple := "ambWT"]
    dt_calls[phased_call == "def_WT", phased_call_simple := "defWT"]

    # Select output columns
    dt_output <- dt_calls[, .(cbc, has_np_data, mutation_call, major_call_simple, phased_call, phased_call_simple)]

    # Add original barcode format (with -1 suffix for 10x)
    dt_output[, barcode_10x := paste0(cbc, "-1")]
    setcolorder(dt_output, c("cbc", "barcode_10x", "has_np_data", "mutation_call", "major_call_simple", "phased_call", "phased_call_simple"))

    # =========================================================================
    # SAVE OUTPUTS
    # =========================================================================

    if (!is.null(path_output_folder)) {
        stopifnot("Output folder not found" = dir.exists(path_output_folder))

        message(Sys.time(), " - Saving outputs")

        out_calls <- file.path(path_output_folder, paste0(session_name, "_flags.csv"))
        fwrite(dt_output, out_calls)
        message("  Saved: ", out_calls)
    }

    return(dt_output)
}
