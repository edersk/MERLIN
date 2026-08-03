# =============================================================================
# summarize_snv.R - Summarize SNV detection results (MERLIN)
# =============================================================================
#
# Takes the per-read table written by detect_snv(), corrects UMIs within each
# cell barcode using Levenshtein (edit) distance, and aggregates reads into a
# per-cell summary. It also detects knee/inflection thresholds on the reads-per-
# UMI distribution and reports SNP allele frequencies so the user can judge
# whether phasing is feasible. Mutation calling itself is left to flag_snv().
# =============================================================================


# --- Internal helpers --------------------------------------------------------

#' Format a count as "n (pct%)" for logging
#' @noRd
.fmt_pct <- function(n, total, digits = 1) {
    paste0(n, " (", round(100 * n / total, digits), "%)")
}


#' Collapse near-identical UMIs within one cell barcode
#'
#' Builds a pairwise Levenshtein distance matrix over the unique UMIs of a cell
#' and greedily lets the most frequent UMI absorb any less-frequent UMI within
#' the edit-distance threshold. Levenshtein distance covers substitutions,
#' insertions and deletions, so it corrects both sequencing errors and indels.
#'
#' @param dt_cbc data.table for a single cell barcode with an \code{initial_UMI}
#'   column.
#' @param umi_mismatch Maximum edit distance for two UMIs to be merged.
#' @param n_cores Unused here; kept for signature compatibility.
#' @return A copy of \code{dt_cbc} with \code{corrected_UMI} and \code{corrected}
#'   columns added.
#' @noRd
.correct_umis_optimized <- function(dt_cbc, umi_mismatch, n_cores) {
    # .SD is locked, so operate on a copy before adding columns by reference.
    dt_cbc <- copy(dt_cbc)

    if (nrow(dt_cbc) <= 1L) {
        dt_cbc[, `:=`(corrected_UMI = initial_UMI, corrected = FALSE)]
        return(dt_cbc)
    }

    # Count and sort UMIs by abundance so the most frequent is processed first.
    umi_counts <- dt_cbc[, .N, by = initial_UMI]
    setorder(umi_counts, -N)

    unique_umis <- umi_counts$initial_UMI
    umi_freq <- umi_counts$N
    n_unique <- length(unique_umis)

    # A single unique UMI needs no correction.
    if (n_unique <= 1L) {
        dt_cbc[, `:=`(corrected_UMI = initial_UMI, corrected = FALSE)]
        return(dt_cbc)
    }

    # Each UMI initially maps to itself.
    umi_to_corrected <- setNames(unique_umis, unique_umis)

    # --- Pairwise Levenshtein distance matrix (with indels) ---
    if (requireNamespace("stringdist", quietly = TRUE)) {
        # method = "lv" is Levenshtein: substitutions + insertions + deletions.
        dist_matrix <- stringdist::stringdistmatrix(unique_umis, unique_umis, method = "lv")
    } else {
        # Fallback edit distance via global alignment if stringdist is missing.
        dist_matrix <- matrix(0L, nrow = n_unique, ncol = n_unique)

        for (i in 1:(n_unique - 1)) {
            for (j in (i + 1):n_unique) {
                aln <- pwalign::pairwiseAlignment(
                    unique_umis[i], unique_umis[j],
                    type = "global",
                    substitutionMatrix = nucleotideSubstitutionMatrix(match = 0, mismatch = 1),
                    gapOpening = 0, gapExtension = 1
                )
                dist <- -pwalign::score(aln) # negative score = edit distance
                dist_matrix[i, j] <- dist
                dist_matrix[j, i] <- dist
            }
        }
    }

    # --- Greedy absorption: frequent UMI swallows nearby rarer UMIs ---
    absorbed <- rep(FALSE, n_unique)

    # Because UMIs are sorted by descending frequency, a lower index means a more
    # frequent UMI; each one absorbs the less-frequent, unabsorbed UMIs near it.
    for (i in seq_len(n_unique)) {
        if (absorbed[i]) next # already merged into a more frequent UMI

        distances_to_i <- dist_matrix[i, ]

        candidates <- which(
            distances_to_i <= umi_mismatch &
                distances_to_i > 0 & # exclude self
                !absorbed &
                seq_len(n_unique) > i # only less-frequent UMIs
        )

        if (length(candidates) > 0) {
            umi_to_corrected[unique_umis[candidates]] <- unique_umis[i]
            absorbed[candidates] <- TRUE
        }
    }

    # Apply the mapping and flag which reads had their UMI changed.
    dt_cbc[, corrected_UMI := umi_to_corrected[initial_UMI]]
    dt_cbc[, corrected := initial_UMI != corrected_UMI]

    return(dt_cbc)
}


#' Detect knee / inflection thresholds on a rank-abundance curve
#'
#' Estimates candidate reads-per-UMI cutoffs using several complementary methods
#' (geometric knee, spline-derivative inflection, lower knee, the \code{uik}
#' method, and a CellRanger-style order-of-magnitude rule), all computed in
#' log-log space. Returns each estimate plus the min/max of the finite ones as a
#' suggested threshold range.
#'
#' @param x Ranks (x-axis of the rank-abundance curve).
#' @param y Reads per UMI (y-axis), aligned to \code{x}.
#' @return Named list of rounded threshold estimates.
#' @noRd
.detect_knee_points <- function(x, y) {
    if (length(x) < 10L) {
        return(list(
            upper_knee = NA_real_, inflection = NA_real_, lower_knee = NA_real_,
            cellranger = NA_real_, uik = NA_real_,
            lower = NA_real_, upper = NA_real_
        ))
    }

    # Keep only strictly positive, finite points (log space requires this).
    valid <- is.finite(x) & is.finite(y) & y > 0 & x > 0
    x <- x[valid]
    y <- y[valid]
    n <- length(x) # total points (used by the CellRanger rule)

    if (n < 10L) {
        return(list(
            upper_knee = NA_real_, inflection = NA_real_, lower_knee = NA_real_,
            cellranger = NA_real_, uik = NA_real_,
            lower = NA_real_, upper = NA_real_
        ))
    }

    # Rank-abundance curves are analysed in log-log space.
    x_log <- log10(x)
    y_log <- log10(y)

    # --- Method 1: geometric knee ---
    # Point of maximum perpendicular distance from the line joining the two
    # endpoints, i.e. where the curve bends most sharply.
    .find_geometric_knee <- function(x_vals, y_vals) {
        n <- length(x_vals)
        if (n < 3) {
            return(NA_integer_)
        }

        x1 <- x_vals[1]
        y1 <- y_vals[1]
        x2 <- x_vals[n]
        y2 <- y_vals[n]

        # Perpendicular distance of each point to the endpoint-to-endpoint line.
        a <- y2 - y1
        b <- -(x2 - x1)
        c <- (x2 - x1) * y1 - (y2 - y1) * x1

        distances <- abs(a * x_vals + b * y_vals + c) / sqrt(a^2 + b^2)

        # Ignore the endpoints themselves.
        distances[1] <- 0
        distances[n] <- 0

        which.max(distances)
    }

    # Only strip obvious noise (UMIs seen 1-2 times); keep enough points to fit.
    meaningful_idx <- which(y > 2)
    n_meaningful <- length(meaningful_idx)

    if (n_meaningful < 20) {
        meaningful_idx <- seq_len(n)
        n_meaningful <- n
    }

    x_meaningful <- x_log[meaningful_idx]
    y_meaningful_raw <- y_log[meaningful_idx]
    y_orig_meaningful <- y[meaningful_idx]

    # --- Smooth the curve (DropletUtils-style) for stable derivatives ---
    spline_fit <- tryCatch(
        {
            smooth.spline(x_meaningful, y_meaningful_raw, df = min(20, n_meaningful / 3))
        },
        error = function(e) NULL
    )

    if (!is.null(spline_fit)) {
        y_fitted <- predict(spline_fit, x_meaningful)$y

        # First derivative (slope) and second derivative (curvature).
        d1 <- predict(spline_fit, x_meaningful, deriv = 1)$y
        d2 <- predict(spline_fit, x_meaningful, deriv = 2)$y
    } else {
        # If the spline fails, fall back to a loess fit and finite differences.
        y_fitted <- tryCatch(
            {
                lo <- loess(y_meaningful_raw ~ x_meaningful, span = 0.1)
                predict(lo)
            },
            error = function(e) y_meaningful_raw
        )

        d1 <- c(0, diff(y_fitted) / diff(x_meaningful))
        d2 <- c(0, diff(d1))
    }

    # --- Knee point (max negative signed curvature) ---
    curvature <- d2 / (1 + d1^2)^1.5

    # Search only the upper half of the curve to avoid the noisy tail.
    upper_portion <- 1:min(floor(n_meaningful * 0.5), n_meaningful)
    knee_idx <- upper_portion[which.min(curvature[upper_portion])]
    upper_knee <- if (length(knee_idx) > 0 && knee_idx >= 1 && knee_idx <= length(y_orig_meaningful)) {
        y_orig_meaningful[knee_idx]
    } else {
        NA_real_
    }

    # --- Inflection point (steepest descent = minimum first derivative) ---
    inflection_idx <- which.min(d1)
    inflection_knee <- if (length(inflection_idx) > 0 && inflection_idx >= 1 && inflection_idx <= length(y_orig_meaningful)) {
        y_orig_meaningful[inflection_idx]
    } else {
        NA_real_
    }

    # --- Lower knee (where the curve flattens: max positive curvature) ---
    lower_portion <- floor(n_meaningful * 0.5):n_meaningful
    lower_knee_local <- lower_portion[which.max(curvature[lower_portion])]
    lower_knee <- if (length(lower_knee_local) > 0 && lower_knee_local >= 1 && lower_knee_local <= length(y_orig_meaningful)) {
        y_orig_meaningful[lower_knee_local]
    } else {
        NA_real_
    }

    # --- UIK (unit-invariant knee) from the inflection package, as a backup ---
    uik_knee <- tryCatch(
        {
            knee_x <- inflection::uik(x_meaningful, y_meaningful_raw)
            knee_idx <- which.min(abs(x_meaningful - knee_x))
            if (knee_idx >= 1 && knee_idx <= length(y_orig_meaningful)) {
                y_orig_meaningful[knee_idx]
            } else {
                NA_real_
            }
        },
        error = function(e) NA_real_
    )

    # --- CellRanger-style order-of-magnitude rule ---
    # 10x uses m = 99th percentile of the top-N barcodes, threshold = m / 10.
    cellranger_knee <- tryCatch(
        {
            n_expected <- max(100, floor(n * 0.01))
            top_n <- min(n_expected, n)
            m <- quantile(y[1:top_n], 0.99)
            round(m / 10)
        },
        error = function(e) NA_real_
    )

    # Combine all finite estimates into a suggested threshold range.
    all_knees <- c(upper_knee, inflection_knee, lower_knee, uik_knee, cellranger_knee)
    all_knees <- all_knees[is.finite(all_knees) & all_knees > 1]

    if (length(all_knees) > 0L) {
        lower <- min(all_knees)
        upper <- max(all_knees)
    } else {
        lower <- upper <- NA_real_
    }

    return(list(
        upper_knee = round(upper_knee),
        inflection = round(inflection_knee),
        lower_knee = round(lower_knee),
        uik = round(uik_knee),
        cellranger = round(cellranger_knee),
        lower = round(lower),
        upper = round(upper)
    ))
}


#' Write a 2x2 diagnostic PDF for the reads-per-UMI distribution
#'
#' Panels: (1) smoothed rank-abundance curve with threshold lines, (2) histogram
#' of reads per UMI, (3) cumulative read distribution, and (4) violin plots of
#' major reads split by minor-ROI classification. Wrapped in tryCatch so a
#' plotting failure never aborts the summarisation.
#'
#' @param umi_freq data.table with \code{rank} and \code{n_reads}.
#' @param knee_points Output of \code{.detect_knee_points}.
#' @param output_folder,session_name Where and how to name the PDF.
#' @param dt_umi_counts Optional per-UMI counts used for the 4th panel.
#' @return Invisibly \code{NULL}; called for the side effect of writing a PDF.
#' @noRd
.plot_umi_distribution <- function(umi_freq, knee_points, output_folder, session_name,
                                   dt_umi_counts = NULL) {
    tryCatch(
        {
            plot_file <- file.path(output_folder, paste0(session_name, "_umi_distribution.pdf"))

            pdf(plot_file, width = 14, height = 10)

            # 2x2 panel layout.
            par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

            # Smooth the rank-abundance data in log-log space for panel 1.
            x_log <- log10(umi_freq$rank)
            y_log <- log10(umi_freq$n_reads)

            smooth_span <- min(0.1, 100 / nrow(umi_freq))
            lo <- tryCatch(
                {
                    loess(y_log ~ x_log, span = smooth_span)
                },
                error = function(e) NULL
            )

            if (!is.null(lo)) {
                y_smooth <- predict(lo)
            } else {
                y_smooth <- y_log # fall back to the raw curve
            }

            # === Panel 1: smoothed reads-per-UMI rank-abundance curve ===
            plot(umi_freq$rank, umi_freq$n_reads,
                type = "n", # empty: draw axes only, then add the smoothed line
                log = "xy",
                xlab = "UMI Rank (log scale)",
                ylab = "Reads per UMI (log scale)",
                main = "Reads per UMI Distribution - Smoothed"
            )

            lines(10^x_log, 10^y_smooth, col = "steelblue", lwd = 2)

            # Overlay each detected threshold as a labelled horizontal line.
            if (is.finite(knee_points$upper_knee)) {
                abline(h = knee_points$upper_knee, col = "red", lwd = 2, lty = 1)
                text(par("usr")[1] + 0.02 * diff(par("usr")[1:2]), knee_points$upper_knee,
                    sprintf("Upper Knee: %d", knee_points$upper_knee),
                    pos = 3, col = "red", cex = 0.8
                )
            }
            if (is.finite(knee_points$inflection)) {
                abline(h = knee_points$inflection, col = "orange", lwd = 2, lty = 2)
                text(par("usr")[1] + 0.02 * diff(par("usr")[1:2]), knee_points$inflection,
                    sprintf("Inflection: %d", knee_points$inflection),
                    pos = 3, col = "orange", cex = 0.8
                )
            }
            if (is.finite(knee_points$lower_knee)) {
                abline(h = knee_points$lower_knee, col = "purple", lwd = 2, lty = 3)
                text(par("usr")[1] + 0.02 * diff(par("usr")[1:2]), knee_points$lower_knee,
                    sprintf("Lower Knee: %d", knee_points$lower_knee),
                    pos = 3, col = "purple", cex = 0.8
                )
            }
            if (is.finite(knee_points$cellranger)) {
                abline(h = knee_points$cellranger, col = "darkgreen", lwd = 2, lty = 4)
                text(par("usr")[1] + 0.02 * diff(par("usr")[1:2]), knee_points$cellranger,
                    sprintf("CellRanger: %d", knee_points$cellranger),
                    pos = 3, col = "darkgreen", cex = 0.8
                )
            }

            legend("topright",
                legend = c(
                    "Reads per UMI (smoothed)",
                    if (is.finite(knee_points$upper_knee)) sprintf("Upper Knee: %d", knee_points$upper_knee) else NULL,
                    if (is.finite(knee_points$inflection)) sprintf("Inflection: %d", knee_points$inflection) else NULL,
                    if (is.finite(knee_points$lower_knee)) sprintf("Lower Knee: %d", knee_points$lower_knee) else NULL,
                    if (is.finite(knee_points$cellranger)) sprintf("CellRanger: %d", knee_points$cellranger) else NULL
                ),
                col = c(
                    "steelblue",
                    if (is.finite(knee_points$upper_knee)) "red" else NULL,
                    if (is.finite(knee_points$inflection)) "orange" else NULL,
                    if (is.finite(knee_points$lower_knee)) "purple" else NULL,
                    if (is.finite(knee_points$cellranger)) "darkgreen" else NULL
                ),
                lwd = 2,
                lty = c(1, 1, 2, 3, 4)[1:(1 + sum(is.finite(c(knee_points$upper_knee, knee_points$inflection, knee_points$lower_knee, knee_points$cellranger))))],
                cex = 0.6, bg = "white"
            )

            # === Panel 2: histogram of reads per UMI (log x) ===
            hist(log10(umi_freq$n_reads + 1),
                breaks = 50,
                col = "lightblue", border = "white",
                xlab = "log10(Reads per UMI)",
                ylab = "Number of UMIs",
                main = "Distribution of Reads per UMI"
            )

            if (is.finite(knee_points$upper_knee)) {
                abline(v = log10(knee_points$upper_knee), col = "red", lwd = 2)
            }
            if (is.finite(knee_points$inflection)) {
                abline(v = log10(knee_points$inflection), col = "orange", lwd = 2, lty = 2)
            }
            if (is.finite(knee_points$cellranger)) {
                abline(v = log10(knee_points$cellranger), col = "darkgreen", lwd = 2, lty = 4)
            }

            # === Panel 3: cumulative read distribution ===
            cumulative <- cumsum(umi_freq$n_reads) / sum(umi_freq$n_reads)
            plot(umi_freq$rank, cumulative,
                type = "l", lwd = 2, col = "darkgreen",
                log = "x",
                xlab = "UMI Rank (log scale)",
                ylab = "Cumulative Fraction of Total Reads",
                main = "Cumulative Read Distribution"
            )

            # Mark where each threshold falls on the cumulative curve.
            if (is.finite(knee_points$upper_knee)) {
                upper_rank <- which.min(abs(umi_freq$n_reads - knee_points$upper_knee))
                abline(v = upper_rank, col = "red", lwd = 2)
                abline(h = cumulative[upper_rank], col = "red", lwd = 1, lty = 3)
            }
            if (is.finite(knee_points$inflection)) {
                infl_rank <- which.min(abs(umi_freq$n_reads - knee_points$inflection))
                abline(v = infl_rank, col = "orange", lwd = 2, lty = 2)
            }
            if (is.finite(knee_points$cellranger)) {
                cr_rank <- which.min(abs(umi_freq$n_reads - knee_points$cellranger))
                abline(v = cr_rank, col = "darkgreen", lwd = 2, lty = 4)
            }

            grid(col = "gray90")

            # === Panel 4: major reads split by minor-ROI classification ===
            if (!is.null(dt_umi_counts) && nrow(dt_umi_counts) > 0) {
                # Keep only UMIs with a reasonable amount of major-ROI support.
                dt_plot <- copy(dt_umi_counts)
                dt_plot[, n_major := n_mut + n_wt]
                dt_plot[, n_minor := n_minor_mut + n_minor_wt]
                dt_plot <- dt_plot[n_major >= 4]

                if (nrow(dt_plot) > 0) {
                    # Major class: whichever allele dominates the UMI.
                    dt_plot[, major_class := ifelse(n_mut > n_wt, "Major MUT", "Major WT")]

                    # Minor class relative to an 80%-of-major support threshold:
                    # "minor NA" if there is too little minor signal to trust.
                    dt_plot[, minor_threshold := n_major * 0.8]
                    dt_plot[, minor_class := ifelse(
                        n_minor == 0 | n_minor < minor_threshold,
                        "minor NA",
                        ifelse(n_minor_mut >= minor_threshold, "minor MUT",
                            ifelse(n_minor_wt >= minor_threshold, "minor WT", "minor NA")
                        )
                    )]

                    # y = log10 of the dominant major count.
                    dt_plot[, y_val := log10(pmax(n_mut, n_wt))]

                    # Group = minor class over major class, in a fixed panel order.
                    dt_plot[, group := paste(minor_class, major_class, sep = "\n")]
                    dt_plot[, group := factor(group, levels = c(
                        "minor NA\nMajor MUT", "minor NA\nMajor WT",
                        "minor MUT\nMajor MUT", "minor MUT\nMajor WT",
                        "minor WT\nMajor MUT", "minor WT\nMajor WT"
                    ))]

                    # Consistent fill colours: MUT red-ish, WT blue-ish.
                    group_colors <- c(
                        "minor NA\nMajor MUT" = "firebrick",
                        "minor NA\nMajor WT" = "steelblue",
                        "minor MUT\nMajor MUT" = "firebrick",
                        "minor MUT\nMajor WT" = "steelblue",
                        "minor WT\nMajor MUT" = "firebrick",
                        "minor WT\nMajor WT" = "steelblue"
                    )

                    group_list <- split(dt_plot$y_val, dt_plot$group)

                    # Drop empty groups but preserve the factor ordering.
                    group_list <- group_list[sapply(group_list, length) > 0]

                    if (length(group_list) > 0) {
                        # Lay out violins with a visual gap between minor groups.
                        n_groups <- length(group_list)
                        group_names <- names(group_list)

                        x_positions <- numeric(n_groups)
                        gap <- 0.5
                        pos <- 1
                        current_minor <- sub("\n.*", "", group_names[1])
                        for (i in seq_along(group_names)) {
                            this_minor <- sub("\n.*", "", group_names[i])
                            if (this_minor != current_minor) {
                                pos <- pos + gap
                                current_minor <- this_minor
                            }
                            x_positions[i] <- pos
                            pos <- pos + 1
                        }

                        y_range <- range(unlist(group_list), na.rm = TRUE)
                        y_range <- c(y_range[1] - 0.1 * diff(y_range), y_range[2] + 0.25 * diff(y_range))

                        plot(NULL,
                            xlim = c(0.5, max(x_positions) + 0.5), ylim = y_range,
                            xlab = "", ylab = "log10(Major Reads)",
                            main = "Major Reads by Minor ROI Classification",
                            xaxt = "n"
                        )

                        # Dashed separators between minor groups.
                        unique_minors <- unique(sub("\n.*", "", group_names))
                        if (length(unique_minors) > 1) {
                            for (i in 2:length(unique_minors)) {
                                idx <- which(sub("\n.*", "", group_names) == unique_minors[i])[1]
                                if (!is.na(idx) && idx > 1) {
                                    sep_x <- (x_positions[idx - 1] + x_positions[idx]) / 2
                                    abline(v = sep_x, col = "gray60", lty = 2, lwd = 1)
                                }
                            }
                        }

                        # Simplified per-violin labels (just MUT / WT).
                        simple_labels <- sub(".*\nMajor ", "", group_names)
                        axis(1, at = x_positions, labels = simple_labels, las = 1, cex.axis = 0.8)

                        # Minor-group headers centred over their violins.
                        for (minor_group in unique_minors) {
                            idx <- which(sub("\n.*", "", group_names) == minor_group)
                            if (length(idx) > 0) {
                                mid_x <- mean(x_positions[idx])
                                label <- gsub("minor ", "", minor_group)
                                text(mid_x, y_range[2] - 0.02 * diff(y_range),
                                    label,
                                    font = 2, cex = 0.9
                                )
                            }
                        }

                        # Draw a density-based violin (or a bar for tiny groups).
                        for (i in seq_along(group_list)) {
                            vals <- group_list[[i]]
                            x_pos <- x_positions[i]
                            if (length(vals) > 2) {
                                dens <- density(vals, bw = "SJ", n = 512)
                                # Scale the density to a fixed half-width of 0.35.
                                dens_scaled <- dens$y / max(dens$y) * 0.35

                                polygon(
                                    x = c(x_pos - dens_scaled, rev(x_pos + dens_scaled)),
                                    y = c(dens$x, rev(dens$x)),
                                    col = group_colors[group_names[i]],
                                    border = "black", lwd = 0.5
                                )
                            } else if (length(vals) > 0) {
                                # Too few points for a density: draw a median tick.
                                segments(x_pos - 0.2, median(vals), x_pos + 0.2, median(vals),
                                    col = group_colors[group_names[i]], lwd = 3
                                )
                            }
                        }

                        # Threshold lines (only meaningful at >= 4 major reads).
                        if (is.finite(knee_points$upper_knee) && knee_points$upper_knee >= 4) {
                            thresh_y <- log10(knee_points$upper_knee)
                            abline(h = thresh_y, col = "red", lwd = 1.5, lty = 1)
                            text(par("usr")[2], thresh_y,
                                sprintf("Upper Knee: %d", knee_points$upper_knee),
                                pos = 2, col = "red", cex = 0.6
                            )
                        }
                        if (is.finite(knee_points$inflection) && knee_points$inflection >= 4) {
                            thresh_y <- log10(knee_points$inflection)
                            abline(h = thresh_y, col = "orange", lwd = 1.5, lty = 2)
                            text(par("usr")[2], thresh_y,
                                sprintf("Inflection: %d", knee_points$inflection),
                                pos = 2, col = "orange", cex = 0.6
                            )
                        }
                        if (is.finite(knee_points$cellranger) && knee_points$cellranger >= 4) {
                            thresh_y <- log10(knee_points$cellranger)
                            abline(h = thresh_y, col = "darkgreen", lwd = 1.5, lty = 4)
                            text(par("usr")[2], thresh_y,
                                sprintf("CellRanger: %d", knee_points$cellranger),
                                pos = 2, col = "darkgreen", cex = 0.6
                            )
                        }

                        # Per-group sample size labels.
                        counts <- sapply(group_list, length)
                        text(seq_along(group_list), par("usr")[4] - 0.1,
                            paste0("n=", counts),
                            cex = 0.6, pos = 1
                        )

                        # Threshold legend.
                        threshold_labels <- c()
                        threshold_colors <- c()
                        threshold_ltys <- c()
                        if (is.finite(knee_points$upper_knee)) {
                            threshold_labels <- c(threshold_labels, sprintf("Upper Knee: %d", knee_points$upper_knee))
                            threshold_colors <- c(threshold_colors, "red")
                            threshold_ltys <- c(threshold_ltys, 1)
                        }
                        if (is.finite(knee_points$inflection)) {
                            threshold_labels <- c(threshold_labels, sprintf("Inflection: %d", knee_points$inflection))
                            threshold_colors <- c(threshold_colors, "orange")
                            threshold_ltys <- c(threshold_ltys, 2)
                        }
                        if (is.finite(knee_points$cellranger)) {
                            threshold_labels <- c(threshold_labels, sprintf("CellRanger: %d", knee_points$cellranger))
                            threshold_colors <- c(threshold_colors, "darkgreen")
                            threshold_ltys <- c(threshold_ltys, 4)
                        }
                        if (length(threshold_labels) > 0) {
                            legend("topright",
                                legend = threshold_labels,
                                col = threshold_colors,
                                lty = threshold_ltys,
                                lwd = 1.5,
                                cex = 0.5, bg = "white",
                                title = "Thresholds"
                            )
                        }
                    } else {
                        plot.new()
                        text(0.5, 0.5, "No UMIs passed filtering", cex = 1.2)
                    }
                } else {
                    plot.new()
                    text(0.5, 0.5, "No UMIs with major count >= 2", cex = 1.2)
                }
            } else {
                # No per-UMI data supplied: leave a placeholder panel.
                plot.new()
                text(0.5, 0.5, "No UMI count data available", cex = 1.2)
            }

            dev.off()

            message("  Saved diagnostic plot: ", plot_file)
        },
        error = function(e) {
            warning("Could not generate plot: ", e$message)
        }
    )
}


#' Build major- and minor-ROI SNP frequency tables and assess phasing
#'
#' Summarises reads/UMIs/cells per allele for the major ROI and, when present,
#' the minor ROI, then decides whether phasing is viable (both minor alleles must
#' be present at reasonable frequency).
#'
#' @param dt_corrected UMI-corrected read table.
#' @param has_minor Whether usable minor-ROI data is present.
#' @return List with \code{major_freq}, \code{minor_freq}, \code{combined},
#'   \code{phasing_viable} and \code{phasing_message}.
#' @noRd
.create_snp_frequency_tables <- function(dt_corrected, has_minor) {
    # --- Major ROI frequency table ---
    # n_cells_with_status counts cells that have >= 1 UMI of that status, so a
    # cell carrying both WT and MUT UMIs is counted under both.
    dt_major_freq <- dt_corrected[, .(
        n_reads = .N,
        n_umis = uniqueN(paste0(CBC, "_", corrected_UMI)),
        n_cells_with_status = uniqueN(CBC)
    ), by = major_status]

    total_reads <- sum(dt_major_freq$n_reads)
    total_umis <- sum(dt_major_freq$n_umis)
    total_cells <- uniqueN(dt_corrected$CBC) # actual unique cells

    dt_major_freq[, `:=`(
        pct_reads = round(100 * n_reads / total_reads, 2),
        pct_umis = round(100 * n_umis / total_umis, 2)
    )]

    setorder(dt_major_freq, -n_reads)
    setnames(dt_major_freq, "major_status", "status")
    dt_major_freq[, roi := "major"]

    # --- Minor ROI frequency table (only if phasing data exists) ---
    dt_minor_freq <- NULL
    phasing_viable <- FALSE
    phasing_message <- "No minor ROI data available"

    if (has_minor) {
        dt_minor <- dt_corrected[!is.na(minor_status) & minor_status %in% c("MUT", "WT")]

        if (nrow(dt_minor) > 0) {
            dt_minor_freq <- dt_minor[, .(
                n_reads = .N,
                n_umis = uniqueN(paste0(CBC, "_", corrected_UMI)),
                n_cells_with_status = uniqueN(CBC)
            ), by = minor_status]

            total_reads_minor <- sum(dt_minor_freq$n_reads)
            total_umis_minor <- sum(dt_minor_freq$n_umis)

            dt_minor_freq[, `:=`(
                pct_reads = round(100 * n_reads / total_reads_minor, 2),
                pct_umis = round(100 * n_umis / total_umis_minor, 2)
            )]

            setorder(dt_minor_freq, -n_reads)
            setnames(dt_minor_freq, "minor_status", "status")
            dt_minor_freq[, roi := "minor"]

            # Phasing needs both minor alleles present at >= 10% of UMIs each.
            has_mut <- "MUT" %in% dt_minor_freq$status
            has_wt <- "WT" %in% dt_minor_freq$status

            if (has_mut && has_wt) {
                pct_mut <- dt_minor_freq[status == "MUT", pct_umis]
                pct_wt <- dt_minor_freq[status == "WT", pct_umis]

                if (min(pct_mut, pct_wt) >= 10) {
                    phasing_viable <- TRUE
                    phasing_message <- sprintf(
                        "PHASING VIABLE: Both alleles detected (MUT: %.1f%%, WT: %.1f%% of UMIs)",
                        pct_mut, pct_wt
                    )
                } else {
                    phasing_message <- sprintf(
                        "PHASING LIMITED: Minor allele underrepresented (MUT: %.1f%%, WT: %.1f%% of UMIs)",
                        pct_mut, pct_wt
                    )
                }
            } else {
                phasing_message <- "PHASING NOT POSSIBLE: Only one allele detected in minor ROI"
            }
        }
    }

    # Stack major + minor for a single combined view.
    if (!is.null(dt_minor_freq)) {
        dt_combined <- rbind(dt_major_freq, dt_minor_freq, fill = TRUE)
    } else {
        dt_combined <- dt_major_freq
    }

    return(list(
        major_freq = dt_major_freq,
        minor_freq = dt_minor_freq,
        combined = dt_combined,
        phasing_viable = phasing_viable,
        phasing_message = phasing_message
    ))
}


# --- Main function -----------------------------------------------------------

#' Summarize SNV detection results (UMI correction, no calling)
#'
#' Processes the per-read table from \code{\link{detect_snv}}, corrects UMIs
#' within each cell barcode using Levenshtein distance, and produces a per-cell
#' summary of reads per corrected UMI. It also detects knee/inflection thresholds
#' on the reads-per-UMI distribution and reports SNP allele frequencies to help
#' decide whether phasing is feasible. Mutation calling is performed separately
#' by \code{\link{flag_snv}}.
#'
#' @param session_name Prefix used for output files.
#' @param path_snv_table Path to the CSV written by \code{\link{detect_snv}}.
#' @param path_barcodes Path to a cell barcode whitelist (\code{.csv} or plain
#'   text); a trailing \code{-1} 10x suffix is handled automatically.
#' @param path_output_folder Existing directory for output files.
#' @param umi_mismatch Maximum edit distance for UMI correction (default: 2).
#' @param n_cores Number of CPU cores (default: 4).
#' @param skip_plots Skip diagnostic plot generation (default: FALSE).
#'
#' @return Invisibly, a list with:
#'   \itemize{
#'     \item \code{summary} -- per-cell UMI counts (input for \code{flag_snv()}).
#'     \item \code{umi_frequencies} -- reads per UMI, ranked, for thresholding.
#'     \item \code{knee_points} -- detected knee/inflection thresholds.
#'     \item \code{snp_frequencies}, \code{phasing_viable}, \code{phasing_message},
#'       \code{has_minor_roi} -- phasing assessment.
#'     \item \code{statistics} -- QC counts.
#'   }
#'
#' @details
#' The per-cell summary CSV is semicolon-separated with columns
#' \code{cbc;umis;n_umi;n_umi_mut;n_umi_wt;n_snp;n_snp_mut;n_snp_wt}, where the
#' count columns are comma-separated per corrected UMI. This is the expected
#' input for \code{\link{flag_snv}}.
#'
#' @seealso \code{\link{detect_snv}}, \code{\link{flag_snv}}
#' @export
summarize_snv <- function(session_name,
                          path_snv_table,
                          path_barcodes,
                          path_output_folder,
                          umi_mismatch = 2L,
                          n_cores = 4L,
                          skip_plots = FALSE) {
    # --- Input validation ---
    stopifnot("`session_name` must be character" = is.character(session_name))
    stopifnot("`path_snv_table` not found" = file.exists(path_snv_table))
    stopifnot("`path_barcodes` not found" = file.exists(path_barcodes))
    stopifnot("`path_output_folder` not found" = dir.exists(path_output_folder))

    # Never request more than (available - 1) cores.
    available_cores <- parallel::detectCores()
    if (n_cores > available_cores - 1L) {
        n_cores <- max(1L, available_cores - 1L)
        message("n_cores adjusted to ", n_cores)
    }

    # --- Load the detect_snv() table ---
    message(Sys.time(), " - Loading SNV table: ", path_snv_table)
    dt_input <- fread(path_snv_table, header = TRUE)
    print(colnames(dt_input))
    n_input <- nrow(dt_input)
    message("  Loaded ", format(n_input, big.mark = ","), " reads")

    req_cols <- c("CBC", "UMI", "major_base", "major_status")
    missing <- setdiff(req_cols, colnames(dt_input))
    if (length(missing) > 0) {
        stop("Missing required columns: ", paste(missing, collapse = ", "))
    }

    # Minor-ROI (phasing) data is optional; only use it if actually populated.
    has_minor <- "minor_status" %in% colnames(dt_input)
    if (has_minor) {
        valid_minor <- !is.na(dt_input$minor_status) &
            dt_input$minor_status %in% c("MUT", "WT")
        if (any(valid_minor)) {
            message("  Found minor ROI data: ", sum(valid_minor), " reads")
        } else {
            has_minor <- FALSE
        }
    }

    # --- Restrict to whitelisted barcodes with a valid major call ---
    message(Sys.time(), " - Loading barcode whitelist")

    # Accept either a .csv or a plain-text list of barcodes.
    if (grepl("\\.csv$", path_barcodes, ignore.case = TRUE)) {
        bc_raw <- fread(path_barcodes, header = FALSE)[[1]]
    } else {
        bc_raw <- readLines(path_barcodes)
    }

    # Use the first 16 bp so a 10x "-1" suffix does not break matching.
    bc_16 <- substr(bc_raw, 1, 16)
    message("  Loaded ", length(bc_16), " barcodes")

    # Keep only reads with a definite MUT/WT major call.
    dt_valid <- dt_input[major_status %in% c("MUT", "WT")]
    n_valid <- nrow(dt_valid)
    message("  Valid major ROI calls: ", .fmt_pct(n_valid, n_input))

    # Keep only reads whose barcode is in the whitelist.
    setkey(dt_valid, CBC)
    dt_valid <- dt_valid[CBC %in% bc_16]
    n_filtered <- nrow(dt_valid)
    message("  After barcode filter: ", .fmt_pct(n_filtered, n_valid))
    message("  Unique CBCs: ", uniqueN(dt_valid$CBC))

    # --- UMI correction, one cell barcode at a time ---
    message(Sys.time(), " - Performing UMI correction (max mismatch = ", umi_mismatch, ")")

    # Trim to the columns we need and rename UMI -> initial_UMI.
    cols_keep <- c("read_id", "CBC", "UMI", "major_status")
    if (has_minor) cols_keep <- c(cols_keep, "minor_status")
    dt_valid <- dt_valid[, ..cols_keep]
    setnames(dt_valid, "UMI", "initial_UMI")

    cols_keep[cols_keep == "UMI"] <- "initial_UMI"

    n_cbcs <- uniqueN(dt_valid$CBC)
    dt_corrected <- dt_valid[,
        {
            n_initial_umis <- uniqueN(initial_UMI)
            result <- .correct_umis_optimized(.SD, umi_mismatch, n_cores)
            n_corrected_umis <- uniqueN(result$corrected_UMI)
            n_collapsed <- n_initial_umis - n_corrected_umis

            if (n_collapsed > 10) { # only log the sizeable corrections
                message(sprintf(
                    "  CBC %d/%d - %d reads | UMIs: %d initial -> %d corrected (%d collapsed, %d%% reduction)",
                    .GRP, n_cbcs, nrow(result),
                    n_initial_umis, n_corrected_umis, n_collapsed,
                    round(100 * n_collapsed / max(n_initial_umis, 1))
                ))
            }

            result
        },
        by = CBC,
        .SDcols = setdiff(cols_keep, "CBC")
    ]

    n_corrected <- sum(dt_corrected$corrected)
    message("  UMIs corrected: ", .fmt_pct(n_corrected, nrow(dt_corrected)))

    rm(dt_valid, dt_input)
    gc(verbose = FALSE)

    # --- Reads-per-UMI distribution and knee detection ---
    message(Sys.time(), " - Computing UMI frequencies (reads per UMI)")

    # Reads per corrected UMI (per cell) = the duplication level.
    umi_freq <- dt_corrected[, .(
        n_reads = .N
    ), by = .(CBC, corrected_UMI)]

    # Rank UMIs from most to least sequenced.
    setorder(umi_freq, -n_reads)
    umi_freq[, rank := .I]

    message("  Total unique UMIs: ", nrow(umi_freq))
    message("  Max reads per UMI: ", max(umi_freq$n_reads))
    message("  Median reads per UMI: ", median(umi_freq$n_reads))

    message(Sys.time(), " - Detecting knee/elbow/inflection points")
    knee_points <- .detect_knee_points(umi_freq$rank, umi_freq$n_reads)

    message("\n  === Knee Point Detection Results (Reads per UMI) ===")
    if (is.finite(knee_points$upper_knee)) message("  Upper Knee (start of steep drop): ", knee_points$upper_knee, " reads")
    if (is.finite(knee_points$inflection)) message("  Inflection Point (geometric): ", knee_points$inflection, " reads")
    if (is.finite(knee_points$lower_knee)) message("  Lower Knee (curve flattens): ", knee_points$lower_knee, " reads")
    if (is.finite(knee_points$uik)) message("  UIK (Unit Invariant Knee): ", knee_points$uik, " reads")
    if (is.finite(knee_points$cellranger)) message("  CellRanger-style (m/10): ", knee_points$cellranger, " reads")
    if (is.finite(knee_points$lower)) {
        message("  Suggested threshold range: ", knee_points$lower, " - ", knee_points$upper, " reads per UMI")
    }
    message("")

    # (Plots are drawn later, once dt_umi_counts exists for the 4th panel.)

    # --- Per-UMI MUT/WT (and minor) counts ---
    message(Sys.time(), " - Summarizing reads per UMI")

    # Count MUT vs WT reads for each corrected UMI.
    dt_umi_counts <- dt_corrected[, .(
        n_mut = sum(major_status == "MUT"),
        n_wt = sum(major_status == "WT"),
        n_total = .N
    ), by = .(CBC, corrected_UMI)]

    # Attach minor-ROI counts when phasing data is available.
    if (has_minor) {
        dt_minor <- dt_corrected[!is.na(minor_status) & minor_status %in% c("MUT", "WT")]

        if (nrow(dt_minor) > 0) {
            dt_minor_counts <- dt_minor[, .(
                n_minor_mut = sum(minor_status == "MUT"),
                n_minor_wt = sum(minor_status == "WT")
            ), by = .(CBC, corrected_UMI)]

            setkey(dt_umi_counts, CBC, corrected_UMI)
            setkey(dt_minor_counts, CBC, corrected_UMI)
            dt_umi_counts <- merge(dt_umi_counts, dt_minor_counts, all.x = TRUE)
            dt_umi_counts[is.na(n_minor_mut), n_minor_mut := 0L]
            dt_umi_counts[is.na(n_minor_wt), n_minor_wt := 0L]
        } else {
            dt_umi_counts[, `:=`(n_minor_mut = 0L, n_minor_wt = 0L)]
        }
    } else {
        dt_umi_counts[, `:=`(n_minor_mut = 0L, n_minor_wt = 0L)]
    }

    # Now that per-UMI counts exist, draw the diagnostic plots.
    if (!skip_plots) {
        .plot_umi_distribution(umi_freq, knee_points, path_output_folder, session_name, dt_umi_counts)
    }

    # Order UMIs within each cell from most to least abundant.
    setorder(dt_umi_counts, CBC, -n_total)

    # Collapse each cell's UMIs into comma-separated strings (one row per cell).
    dt_summary <- dt_umi_counts[, .(
        umis = paste(corrected_UMI, collapse = ","),
        n_umi = paste(n_total, collapse = ","),
        n_umi_mut = paste(n_mut, collapse = ","),
        n_umi_wt = paste(n_wt, collapse = ","),
        n_snp = paste(n_minor_mut + n_minor_wt, collapse = ","),
        n_snp_mut = paste(n_minor_mut, collapse = ","),
        n_snp_wt = paste(n_minor_wt, collapse = ",")
    ), by = CBC]

    setnames(dt_summary, "CBC", "cbc")

    message("  Summary created for ", nrow(dt_summary), " cell barcodes")

    # --- Write outputs ---
    message(Sys.time(), " - Saving outputs")

    # Per-cell summary (semicolon-separated) = input for flag_snv().
    out_summary <- file.path(path_output_folder, paste0(session_name, "_snv_summary.csv"))
    fwrite(dt_summary, out_summary, sep = ";")
    message("  ", out_summary)

    # Ranked reads-per-UMI table (for choosing a threshold).
    out_freq <- file.path(path_output_folder, paste0(session_name, "_umi_frequencies.csv"))
    fwrite(umi_freq, out_freq)
    message("  ", out_freq)

    # --- SNP frequency tables + phasing assessment ---
    message(Sys.time(), " - Creating SNP frequency tables")

    snp_freq_tables <- .create_snp_frequency_tables(dt_corrected, has_minor)

    out_major_freq <- file.path(path_output_folder, paste0(session_name, "_major_roi_frequency.csv"))
    fwrite(snp_freq_tables$major_freq, out_major_freq)
    message("  ", out_major_freq)

    if (!is.null(snp_freq_tables$minor_freq)) {
        out_minor_freq <- file.path(path_output_folder, paste0(session_name, "_minor_roi_frequency.csv"))
        fwrite(snp_freq_tables$minor_freq, out_minor_freq)
        message("  ", out_minor_freq)
    }

    out_combined_freq <- file.path(path_output_folder, paste0(session_name, "_snp_frequency_combined.csv"))
    fwrite(snp_freq_tables$combined, out_combined_freq)
    message("  ", out_combined_freq)

    # Report the phasing assessment to the console.
    message("\n  === SNP Frequency Summary ===")
    message("  Major ROI:")
    print(snp_freq_tables$major_freq[, .(status, n_reads, pct_reads, n_umis, pct_umis)])

    if (!is.null(snp_freq_tables$minor_freq)) {
        message("\n  Minor ROI:")
        print(snp_freq_tables$minor_freq[, .(status, n_reads, pct_reads, n_umis, pct_umis)])
    }

    message("\n  ", snp_freq_tables$phasing_message)

    # --- Return ---
    message("\n=== Summarization Complete ===")
    message("Input reads: ", format(n_input, big.mark = ","))
    message("Valid reads: ", format(n_filtered, big.mark = ","))
    message("Cells summarized: ", nrow(dt_summary))
    message("UMIs corrected: ", n_corrected)
    message("\nNext step: Use flag_snv() to perform mutation calling")

    result <- list(
        summary = dt_summary,
        umi_frequencies = umi_freq,
        knee_points = knee_points,
        snp_frequencies = snp_freq_tables,
        phasing_viable = snp_freq_tables$phasing_viable,
        phasing_message = snp_freq_tables$phasing_message,
        has_minor_roi = has_minor,
        statistics = list(
            n_input = n_input,
            n_valid = n_valid,
            n_filtered = n_filtered,
            n_cells = nrow(dt_summary),
            n_umis_corrected = n_corrected
        )
    )

    invisible(result)
}
