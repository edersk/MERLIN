# =============================================================================
# summarize_snv.R - SNV Detection Results Summarization
# =============================================================================
#
# Description:
#   Processes output from detect_snv(), performs UMI error correction,
#   and creates per-cell barcode summaries. Does NOT perform mutation calling.
#   Use flag_snv() for mutation calling with thresholds.
#
# Author: Sebastian Eder <sebastian.eder@example.com>
# Created: 2025-12-18
# Version: 2.1.0
#
# Pipeline:
#   detect_snv() -> summarize_snv() -> flag_snv()
#
# Features:
#   - UMI correction for insertions, deletions, and substitutions
#   - Per-CBC UMI summary with read counts
#   - Knee/elbow/inflection point detection and diagnostic plots
#   - Output compatible with flag_snv() for mutation calling
#
# License: MIT License
#   Copyright (c) 2025 Sebastian Eder
#
# Dependencies:
#   - Biostrings (>= 2.60.0)
#   - data.table (>= 1.14.0)
#   - parallel (>= 4.0.0)
#   - inflection (>= 1.3.5)
#
# =============================================================================


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Formats percentage for logging
.fmt_pct <- function(n, total, digits = 1) {
    paste0(n, " (", round(100 * n / total, digits), "%)")
}


# Corrects UMI sequences within each CBC using pairwise distance matrix
# Uses Levenshtein distance (edit distance) including insertions, deletions, substitutions
# Algorithm: Most frequent UMI "absorbs" less frequent UMIs within distance threshold
.correct_umis_optimized <- function(dt_cbc, umi_mismatch, n_cores) {
    # Make a copy - .SD is locked and cannot be modified with :=
    dt_cbc <- copy(dt_cbc)

    if (nrow(dt_cbc) <= 1L) {
        dt_cbc[, `:=`(corrected_UMI = initial_UMI, corrected = FALSE)]
        return(dt_cbc)
    }

    # Count UMI frequencies and sort by abundance (most frequent first)
    umi_counts <- dt_cbc[, .N, by = initial_UMI]
    setorder(umi_counts, -N)

    unique_umis <- umi_counts$initial_UMI
    umi_freq <- umi_counts$N
    n_unique <- length(unique_umis)

    # If only 1 unique UMI, nothing to correct
    if (n_unique <= 1L) {
        dt_cbc[, `:=`(corrected_UMI = initial_UMI, corrected = FALSE)]
        return(dt_cbc)
    }

    # Initialize mapping: each UMI maps to itself initially
    umi_to_corrected <- setNames(unique_umis, unique_umis)

    # =========================================================================
    # GREEDY UMI CORRECTION (memory-efficient for large CBCs)
    # =========================================================================
    # For small CBCs (< 5000 unique UMIs): full pairwise distance matrix
    # For large CBCs: compare each low-freq UMI against high-freq ones only

    absorbed <- rep(FALSE, n_unique)
    max_matrix_size <- 5000L  # threshold for switching to row-wise computation

    if (n_unique <= max_matrix_size) {
        # Full distance matrix approach (fast for small CBCs)
        if (requireNamespace("stringdist", quietly = TRUE)) {
            dist_matrix <- stringdist::stringdistmatrix(unique_umis, unique_umis, method = "lv")
        } else {
            dist_matrix <- matrix(0L, nrow = n_unique, ncol = n_unique)
            for (i in 1:(n_unique - 1)) {
                for (j in (i + 1):n_unique) {
                    aln <- pwalign::pairwiseAlignment(
                        unique_umis[i], unique_umis[j], type = "global",
                        substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 0, mismatch = 1),
                        gapOpening = 0, gapExtension = 1)
                    dist <- -pwalign::score(aln)
                    dist_matrix[i, j] <- dist
                    dist_matrix[j, i] <- dist
                }
            }
        }

        for (i in seq_len(n_unique)) {
            if (absorbed[i]) next
            distances_to_i <- dist_matrix[i, ]
            candidates <- which(
                distances_to_i <= umi_mismatch & distances_to_i > 0 &
                !absorbed & seq_len(n_unique) > i)
            if (length(candidates) > 0) {
                umi_to_corrected[unique_umis[candidates]] <- unique_umis[i]
                absorbed[candidates] <- TRUE
            }
        }
        rm(dist_matrix)
    } else {
        # Row-wise approach for large CBCs: compute distances one row at a time
        # Process most frequent first; only compare against not-yet-absorbed UMIs
        has_stringdist <- requireNamespace("stringdist", quietly = TRUE)

        for (i in seq_len(n_unique)) {
            if (absorbed[i]) next

            # Only compare against less-frequent, non-absorbed UMIs
            remaining <- which(!absorbed & seq_len(n_unique) > i)
            if (length(remaining) == 0L) next

            if (has_stringdist) {
                dists <- stringdist::stringdist(unique_umis[i], unique_umis[remaining], method = "lv")
            } else {
                dists <- vapply(unique_umis[remaining], function(u) {
                    aln <- pwalign::pairwiseAlignment(unique_umis[i], u, type = "global",
                        substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix(match = 0, mismatch = 1),
                        gapOpening = 0, gapExtension = 1)
                    -pwalign::score(aln)
                }, numeric(1))
            }

            close_idx <- remaining[dists <= umi_mismatch & dists > 0]
            if (length(close_idx) > 0) {
                umi_to_corrected[unique_umis[close_idx]] <- unique_umis[i]
                absorbed[close_idx] <- TRUE
            }
        }
    }

    # Apply corrections
    dt_cbc[, corrected_UMI := umi_to_corrected[initial_UMI]]
    dt_cbc[, corrected := initial_UMI != corrected_UMI]

    return(dt_cbc)
}


# All-NA knee point result (too few points to estimate anything).
# Must contain every key that .detect_knee_points() returns on success so
# downstream is.finite() checks never see NULL.
.empty_knee_points <- function() {
    list(
        upper_knee = NA_real_, inflection = NA_real_, lower_knee = NA_real_,
        uik = NA_real_, cellranger = NA_real_, otsu = NA_real_, mixture = NA_real_,
        lower = NA_real_, upper = NA_real_
    )
}


# Detect knee/elbow points using multiple methods
# Works in log-log space for rank-abundance distributions
.detect_knee_points <- function(x, y) {
    if (length(x) < 10L) {
        return(.empty_knee_points())
    }

    # Remove invalid values
    valid <- is.finite(x) & is.finite(y) & y > 0 & x > 0
    x <- x[valid]
    y <- y[valid]
    n <- length(x) # Total number of points (for CellRanger calculation)

    if (n < 10L) {
        return(.empty_knee_points())
    }

    # Work in log-log space (appropriate for rank-abundance curves)
    x_log <- log10(x)
    y_log <- log10(y)

    # === Method 1: Geometric knee detection ===
    # Find point of maximum perpendicular distance from line connecting endpoints
    # This identifies where the curve bends most sharply

    .find_geometric_knee <- function(x_vals, y_vals) {
        n <- length(x_vals)
        if (n < 3) {
            return(NA_integer_)
        }

        # Line from first to last point
        x1 <- x_vals[1]
        y1 <- y_vals[1]
        x2 <- x_vals[n]
        y2 <- y_vals[n]

        # Calculate perpendicular distance for each point
        # Distance = |ax + by + c| / sqrt(a^2 + b^2)
        # where line is: (y2-y1)x - (x2-x1)y + (x2-x1)y1 - (y2-y1)x1 = 0
        a <- y2 - y1
        b <- -(x2 - x1)
        c <- (x2 - x1) * y1 - (y2 - y1) * x1

        distances <- abs(a * x_vals + b * y_vals + c) / sqrt(a^2 + b^2)

        # Find point with maximum distance (excluding endpoints)
        distances[1] <- 0
        distances[n] <- 0

        which.max(distances)
    }

    # Only exclude obvious noise (y <= 2)
    meaningful_idx <- which(y > 2)
    n_meaningful <- length(meaningful_idx)

    if (n_meaningful < 20) {
        meaningful_idx <- seq_len(n)
        n_meaningful <- n
    }

    x_meaningful <- x_log[meaningful_idx]
    y_meaningful_raw <- y_log[meaningful_idx]
    y_orig_meaningful <- y[meaningful_idx]

    # === Apply smooth.spline for curve fitting (DropletUtils approach) ===
    # Use smooth.spline with controlled degrees of freedom for stable derivatives
    spline_fit <- tryCatch(
        {
            smooth.spline(x_meaningful, y_meaningful_raw, df = min(20, n_meaningful / 3))
        },
        error = function(e) NULL
    )

    if (!is.null(spline_fit)) {
        y_fitted <- predict(spline_fit, x_meaningful)$y

        # Calculate derivatives from spline
        # First derivative (slope)
        d1 <- predict(spline_fit, x_meaningful, deriv = 1)$y

        # Second derivative (curvature)
        d2 <- predict(spline_fit, x_meaningful, deriv = 2)$y
    } else {
        # Fallback: use loess if spline fails
        y_fitted <- tryCatch(
            {
                lo <- loess(y_meaningful_raw ~ x_meaningful, span = 0.1)
                predict(lo)
            },
            error = function(e) y_meaningful_raw
        )

        # Approximate derivatives from smoothed data
        d1 <- c(0, diff(y_fitted) / diff(x_meaningful))
        d2 <- c(0, diff(d1))
    }

    # === KNEE POINT (DropletUtils method) ===
    # Knee = point where signed curvature is minimized (maximum negative curvature)
    # Curvature = d2y / (1 + d1y^2)^(3/2)
    curvature <- d2 / (1 + d1^2)^1.5

    # Find knee: minimum curvature in the upper portion of curve (exclude noisy tail)
    upper_portion <- 1:min(floor(n_meaningful * 0.5), n_meaningful)
    knee_idx <- upper_portion[which.min(curvature[upper_portion])]
    upper_knee <- if (length(knee_idx) > 0 && knee_idx >= 1 && knee_idx <= length(y_orig_meaningful)) {
        y_orig_meaningful[knee_idx]
    } else {
        NA_real_
    }

    # === INFLECTION POINT (DropletUtils method) ===
    # Inflection = point where first derivative is minimized (steepest descent)
    # Look across the full meaningful range
    inflection_idx <- which.min(d1)
    inflection_knee <- if (length(inflection_idx) > 0 && inflection_idx >= 1 && inflection_idx <= length(y_orig_meaningful)) {
        y_orig_meaningful[inflection_idx]
    } else {
        NA_real_
    }

    # === LOWER KNEE (transition to flat) ===
    # Find where curvature becomes most positive (curve flattening) in lower portion
    lower_portion <- floor(n_meaningful * 0.5):n_meaningful
    lower_knee_local <- lower_portion[which.max(curvature[lower_portion])]
    lower_knee <- if (length(lower_knee_local) > 0 && lower_knee_local >= 1 && lower_knee_local <= length(y_orig_meaningful)) {
        y_orig_meaningful[lower_knee_local]
    } else {
        NA_real_
    }

    # === UIK from inflection package (backup method) ===
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

    # === CellRanger-style OrdMag ===
    # 10x Genomics uses: m = 99th percentile of top N barcodes, threshold = m/10
    cellranger_knee <- tryCatch(
        {
            n_expected <- max(100, floor(n * 0.01))
            top_n <- min(n_expected, n)
            m <- quantile(y[1:top_n], 0.99)
            round(m / 10)
        },
        error = function(e) NA_real_
    )

    # === Otsu's method ===
    # Maximize inter-class variance on log10(reads per UMI) histogram
    # Borrowed from image thresholding — finds the cut that best separates
    # two populations (signal vs noise) in the intensity histogram
    otsu_knee <- tryCatch(
        {
            log_y <- log10(y[y > 0])
            n_bins <- min(256L, length(unique(log_y)))
            h <- hist(log_y, breaks = n_bins, plot = FALSE)
            counts <- h$counts
            mids <- h$mids
            total <- sum(counts)

            best_variance <- -Inf
            best_threshold <- NA_real_

            cum_sum_w <- 0
            cum_sum_wm <- 0
            global_mean <- sum(counts * mids) / total

            for (t in seq_len(length(counts) - 1L)) {
                cum_sum_w <- cum_sum_w + counts[t]
                cum_sum_wm <- cum_sum_wm + counts[t] * mids[t]

                w0 <- cum_sum_w / total
                w1 <- 1 - w0
                if (w0 == 0 || w1 == 0) next

                mu0 <- cum_sum_wm / cum_sum_w
                mu1 <- (global_mean * total - cum_sum_wm) / (total - cum_sum_w)

                between_var <- w0 * w1 * (mu0 - mu1)^2

                if (between_var > best_variance) {
                    best_variance <- between_var
                    best_threshold <- mids[t]
                }
            }

            if (!is.na(best_threshold)) round(10^best_threshold) else NA_real_
        },
        error = function(e) NA_real_
    )

    # === Mixture model (2-component Gaussian in log-space) ===
    # Fits two normal distributions to log10(reads per UMI) using EM algorithm
    # The intersection of the two components defines the threshold
    mixture_knee <- tryCatch(
        {
            log_y <- log10(y[y > 0])
            n_obs <- length(log_y)
            if (n_obs < 50L) stop("too few observations")

            # Initialize: split at median
            med <- median(log_y)
            mu1 <- mean(log_y[log_y <= med])
            mu2 <- mean(log_y[log_y > med])
            sd1 <- sd(log_y[log_y <= med])
            sd2 <- sd(log_y[log_y > med])
            pi1 <- 0.5

            if (is.na(sd1) || sd1 == 0) sd1 <- sd(log_y) / 2
            if (is.na(sd2) || sd2 == 0) sd2 <- sd(log_y) / 2

            # EM iterations
            for (iter in seq_len(100L)) {
                # E-step: posterior probability of component 1
                d1 <- pi1 * dnorm(log_y, mu1, sd1)
                d2 <- (1 - pi1) * dnorm(log_y, mu2, sd2)
                total_d <- d1 + d2
                total_d[total_d == 0] <- .Machine$double.xmin
                gamma <- d1 / total_d

                # M-step
                n1 <- sum(gamma)
                n2 <- n_obs - n1
                if (n1 < 2 || n2 < 2) break

                pi1_new <- n1 / n_obs
                mu1_new <- sum(gamma * log_y) / n1
                mu2_new <- sum((1 - gamma) * log_y) / n2
                sd1_new <- sqrt(sum(gamma * (log_y - mu1_new)^2) / n1)
                sd2_new <- sqrt(sum((1 - gamma) * (log_y - mu2_new)^2) / n2)

                if (sd1_new < 1e-6) sd1_new <- 1e-6
                if (sd2_new < 1e-6) sd2_new <- 1e-6

                # Check convergence
                if (abs(mu1_new - mu1) + abs(mu2_new - mu2) < 1e-6) break

                mu1 <- mu1_new; mu2 <- mu2_new
                sd1 <- sd1_new; sd2 <- sd2_new
                pi1 <- pi1_new
            }

            # Ensure mu1 < mu2 (component 1 = noise, component 2 = signal)
            if (mu1 > mu2) {
                tmp <- mu1; mu1 <- mu2; mu2 <- tmp
                tmp <- sd1; sd1 <- sd2; sd2 <- tmp
                pi1 <- 1 - pi1
            }

            # Find intersection between mu1 and mu2
            search_grid <- seq(mu1, mu2, length.out = 1000L)
            f1 <- pi1 * dnorm(search_grid, mu1, sd1)
            f2 <- (1 - pi1) * dnorm(search_grid, mu2, sd2)
            cross_idx <- which(diff(sign(f2 - f1)) != 0)

            if (length(cross_idx) > 0) {
                threshold_log <- search_grid[cross_idx[1]]
                round(10^threshold_log)
            } else {
                # Fallback: midpoint between means
                round(10^((mu1 + mu2) / 2))
            }
        },
        error = function(e) NA_real_
    )

    # Compile results
    all_knees <- c(upper_knee, inflection_knee, lower_knee, uik_knee,
                   cellranger_knee, otsu_knee, mixture_knee)
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
        otsu = round(otsu_knee),
        mixture = round(mixture_knee),
        lower = round(lower),
        upper = round(upper)
    ))
}


# Generate comprehensive diagnostic plots for reads per UMI
.plot_umi_distribution <- function(umi_freq, knee_points, output_folder, session_name,
                                   dt_umi_counts = NULL) {
    tryCatch(
        {
            plot_file <- file.path(output_folder, paste0(session_name, "_umi_distribution.pdf"))

            pdf(plot_file, width = 14, height = 10)

            # Layout: 2x2 grid
            par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

            # Prepare smoothed data for plotting (in log-log space)
            x_log <- log10(umi_freq$rank)
            y_log <- log10(umi_freq$n_reads)

            # Use loess smoothing
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
                y_smooth <- y_log # Fall back to raw data
            }

            # === Plot 1: Ranked reads per UMI distribution (log-log scale) - SMOOTHED ===
            plot(umi_freq$rank, umi_freq$n_reads,
                type = "n", # Empty plot for axes
                log = "xy",
                xlab = "UMI Rank (log scale)",
                ylab = "Reads per UMI (log scale)",
                main = "Reads per UMI Distribution - Smoothed"
            )

            # Add smoothed line
            lines(10^x_log, 10^y_smooth, col = "steelblue", lwd = 2)

            # Add knee point lines with labels
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
            if (is.finite(knee_points$otsu)) {
                abline(h = knee_points$otsu, col = "brown", lwd = 2, lty = 5)
                text(par("usr")[1] + 0.02 * diff(par("usr")[1:2]), knee_points$otsu,
                    sprintf("Otsu: %d", knee_points$otsu),
                    pos = 3, col = "brown", cex = 0.8
                )
            }
            if (is.finite(knee_points$mixture)) {
                abline(h = knee_points$mixture, col = "deeppink", lwd = 2, lty = 6)
                text(par("usr")[1] + 0.02 * diff(par("usr")[1:2]), knee_points$mixture,
                    sprintf("Mixture: %d", knee_points$mixture),
                    pos = 3, col = "deeppink", cex = 0.8
                )
            }

            # Build legend dynamically
            leg_labels <- "Reads per UMI (smoothed)"
            leg_cols <- "steelblue"
            leg_ltys <- 1L
            kp_entries <- list(
                list(knee_points$upper_knee, "Upper Knee", "red", 1L),
                list(knee_points$inflection, "Inflection", "orange", 2L),
                list(knee_points$lower_knee, "Lower Knee", "purple", 3L),
                list(knee_points$cellranger, "CellRanger", "darkgreen", 4L),
                list(knee_points$otsu, "Otsu", "brown", 5L),
                list(knee_points$mixture, "Mixture", "deeppink", 6L)
            )
            for (kp in kp_entries) {
                if (is.finite(kp[[1]])) {
                    leg_labels <- c(leg_labels, sprintf("%s: %d", kp[[2]], kp[[1]]))
                    leg_cols <- c(leg_cols, kp[[3]])
                    leg_ltys <- c(leg_ltys, kp[[4]])
                }
            }
            legend("topright",
                legend = leg_labels, col = leg_cols, lwd = 2, lty = leg_ltys,
                cex = 0.6, bg = "white"
            )

            # === Plot 2: Histogram of reads per UMI (log x-axis) ===
            hist(log10(umi_freq$n_reads + 1),
                breaks = 50,
                col = "lightblue", border = "white",
                xlab = "log10(Reads per UMI)",
                ylab = "Number of UMIs",
                main = "Distribution of Reads per UMI"
            )

            # Add knee lines
            if (is.finite(knee_points$upper_knee)) {
                abline(v = log10(knee_points$upper_knee), col = "red", lwd = 2)
            }
            if (is.finite(knee_points$inflection)) {
                abline(v = log10(knee_points$inflection), col = "orange", lwd = 2, lty = 2)
            }
            if (is.finite(knee_points$cellranger)) {
                abline(v = log10(knee_points$cellranger), col = "darkgreen", lwd = 2, lty = 4)
            }
            if (is.finite(knee_points$otsu)) {
                abline(v = log10(knee_points$otsu), col = "brown", lwd = 2, lty = 5)
            }
            if (is.finite(knee_points$mixture)) {
                abline(v = log10(knee_points$mixture), col = "deeppink", lwd = 2, lty = 6)
            }

            # === Plot 3: Cumulative distribution (log x-axis) ===
            cumulative <- cumsum(umi_freq$n_reads) / sum(umi_freq$n_reads)
            plot(umi_freq$rank, cumulative,
                type = "l", lwd = 2, col = "darkgreen",
                log = "x",
                xlab = "UMI Rank (log scale)",
                ylab = "Cumulative Fraction of Total Reads",
                main = "Cumulative Read Distribution"
            )

            # Mark knee points on cumulative
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
            if (is.finite(knee_points$otsu)) {
                otsu_rank <- which.min(abs(umi_freq$n_reads - knee_points$otsu))
                abline(v = otsu_rank, col = "brown", lwd = 2, lty = 5)
            }
            if (is.finite(knee_points$mixture)) {
                mix_rank <- which.min(abs(umi_freq$n_reads - knee_points$mixture))
                abline(v = mix_rank, col = "deeppink", lwd = 2, lty = 6)
            }

            grid(col = "gray90")

            # === Plot 4: Violin plots of major reads by minor ROI classification ===
            if (!is.null(dt_umi_counts) && nrow(dt_umi_counts) > 0) {
                # Filter: only UMIs with major overall count >= 4
                dt_plot <- copy(dt_umi_counts)
                dt_plot[, n_major := n_mut + n_wt]
                dt_plot[, n_minor := n_minor_mut + n_minor_wt]
                dt_plot <- dt_plot[n_major >= 4]

                if (nrow(dt_plot) > 0) {
                    # Classify major status: MUT if n_mut > n_wt, else WT
                    dt_plot[, major_class := ifelse(n_mut > n_wt, "Major MUT", "Major WT")]

                    # Classify minor status based on 80% threshold of major count
                    # minor NA: no minor counts OR minor < 80% of major
                    # minor MUT: minor_mut >= 80% of major
                    # minor WT: minor_wt >= 80% of major
                    dt_plot[, minor_threshold := n_major * 0.8]
                    dt_plot[, minor_class := ifelse(
                        n_minor == 0 | n_minor < minor_threshold,
                        "minor NA",
                        ifelse(n_minor_mut >= minor_threshold, "minor MUT",
                            ifelse(n_minor_wt >= minor_threshold, "minor WT", "minor NA")
                        )
                    )]

                    # Y values: log10 of major reads (use dominant count)
                    dt_plot[, y_val := log10(pmax(n_mut, n_wt))]

                    # Create factor for grouping - ordered by minor status first
                    dt_plot[, group := paste(minor_class, major_class, sep = "\n")]
                    dt_plot[, group := factor(group, levels = c(
                        "minor NA\nMajor MUT", "minor NA\nMajor WT",
                        "minor MUT\nMajor MUT", "minor MUT\nMajor WT",
                        "minor WT\nMajor MUT", "minor WT\nMajor WT"
                    ))]

                    # Colors for groups - consistent colors for Major MUT vs WT
                    group_colors <- c(
                        "minor NA\nMajor MUT" = "firebrick",
                        "minor NA\nMajor WT" = "steelblue",
                        "minor MUT\nMajor MUT" = "firebrick",
                        "minor MUT\nMajor WT" = "steelblue",
                        "minor WT\nMajor MUT" = "firebrick",
                        "minor WT\nMajor WT" = "steelblue"
                    )

                    # Prepare data for violin plots — keep ALL defined groups
                    all_group_levels <- levels(dt_plot$group)
                    group_list <- split(dt_plot$y_val, dt_plot$group)
                    # Ensure all levels exist (empty ones get length-0 vectors)
                    for (gl in all_group_levels) {
                        if (is.null(group_list[[gl]])) group_list[[gl]] <- numeric(0)
                    }
                    group_list <- group_list[all_group_levels]

                    # Check there is any data at all
                    has_any_data <- any(sapply(group_list, length) > 0)

                    if (has_any_data) {
                        n_groups <- length(group_list)
                        group_names <- names(group_list)

                        # Calculate x positions with gaps between minor groups
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

                        non_empty_vals <- unlist(group_list[sapply(group_list, length) > 0])
                        y_range <- range(non_empty_vals, na.rm = TRUE)
                        y_range <- c(y_range[1] - 0.15 * diff(y_range), y_range[2] + 0.25 * diff(y_range))

                        # Set up empty plot
                        plot(NULL,
                            xlim = c(0.5, max(x_positions) + 0.5), ylim = y_range,
                            xlab = "", ylab = "log10(Major Reads)",
                            main = "Major Reads by Minor ROI Classification",
                            xaxt = "n"
                        )

                        # Add group separators
                        unique_minors <- unique(sub("\n.*", "", group_names))
                        if (length(unique_minors) > 1) {
                            for (i in 2:length(unique_minors)) {
                                # Find where the gap is
                                idx <- which(sub("\n.*", "", group_names) == unique_minors[i])[1]
                                if (!is.na(idx) && idx > 1) {
                                    sep_x <- (x_positions[idx - 1] + x_positions[idx]) / 2
                                    abline(v = sep_x, col = "gray60", lty = 2, lwd = 1)
                                }
                            }
                        }

                        # Add x-axis labels (simplified: just MUT/WT)
                        simple_labels <- sub(".*\nMajor ", "", group_names)
                        axis(1, at = x_positions, labels = simple_labels, las = 1, cex.axis = 0.8)

                        # Add minor group labels at top
                        for (minor_group in unique_minors) {
                            idx <- which(sub("\n.*", "", group_names) == minor_group)
                            if (length(idx) > 0) {
                                mid_x <- mean(x_positions[idx])
                                # Display minor group label
                                label <- gsub("minor ", "", minor_group)
                                text(mid_x, y_range[2] - 0.02 * diff(y_range),
                                    label,
                                    font = 2, cex = 0.9
                                )
                            }
                        }

                        # Draw violin for each group
                        for (i in seq_along(group_list)) {
                            vals <- group_list[[i]]
                            x_pos <- x_positions[i]
                            if (length(vals) > 2) {
                                # Compute density
                                dens <- density(vals, bw = "SJ", n = 512)
                                # Scale density to fit within 0.35 width
                                dens_scaled <- dens$y / max(dens$y) * 0.35

                                # Draw polygon (violin shape)
                                polygon(
                                    x = c(x_pos - dens_scaled, rev(x_pos + dens_scaled)),
                                    y = c(dens$x, rev(dens$x)),
                                    col = group_colors[group_names[i]],
                                    border = "black", lwd = 0.5
                                )
                            } else if (length(vals) > 0) {
                                # For small n, just draw a horizontal line
                                segments(x_pos - 0.2, median(vals), x_pos + 0.2, median(vals),
                                    col = group_colors[group_names[i]], lwd = 3
                                )
                            }
                        }

                        # Add threshold lines with labels
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
                        if (is.finite(knee_points$otsu) && knee_points$otsu >= 4) {
                            thresh_y <- log10(knee_points$otsu)
                            abline(h = thresh_y, col = "brown", lwd = 1.5, lty = 5)
                            text(par("usr")[2], thresh_y,
                                sprintf("Otsu: %d", knee_points$otsu),
                                pos = 2, col = "brown", cex = 0.6
                            )
                        }
                        if (is.finite(knee_points$mixture) && knee_points$mixture >= 4) {
                            thresh_y <- log10(knee_points$mixture)
                            abline(h = thresh_y, col = "deeppink", lwd = 1.5, lty = 6)
                            text(par("usr")[2], thresh_y,
                                sprintf("Mixture: %d", knee_points$mixture),
                                pos = 2, col = "deeppink", cex = 0.6
                            )
                        }

                        # Add count labels at the bottom of each group
                        counts <- sapply(group_list, length)
                        text(x_positions, y_range[1] + 0.02 * diff(y_range),
                            paste0("n=", counts),
                            cex = 0.6, pos = 3
                        )

                        # Add threshold legend
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
                        if (is.finite(knee_points$otsu)) {
                            threshold_labels <- c(threshold_labels, sprintf("Otsu: %d", knee_points$otsu))
                            threshold_colors <- c(threshold_colors, "brown")
                            threshold_ltys <- c(threshold_ltys, 5)
                        }
                        if (is.finite(knee_points$mixture)) {
                            threshold_labels <- c(threshold_labels, sprintf("Mixture: %d", knee_points$mixture))
                            threshold_colors <- c(threshold_colors, "deeppink")
                            threshold_ltys <- c(threshold_ltys, 6)
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
                # Placeholder if no data
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


# Create SNP frequency tables for major and minor ROI
# Helps user decide if phasing is possible (minor ROI needs both alleles)
.create_snp_frequency_tables <- function(dt_corrected, has_minor) {
    # === Major ROI frequency table ===
    # Note: n_cells_with_status counts cells having UMIs of that status
    # A cell with both WT and MUT UMIs will appear in BOTH counts
    dt_major_freq <- dt_corrected[, .(
        n_reads = .N,
        n_umis = uniqueN(paste0(CBC, "_", corrected_UMI)),
        n_cells_with_status = uniqueN(CBC)
    ), by = major_status]

    # Add percentages
    total_reads <- sum(dt_major_freq$n_reads)
    total_umis <- sum(dt_major_freq$n_umis)
    total_cells <- uniqueN(dt_corrected$CBC) # Actual unique cells

    dt_major_freq[, `:=`(
        pct_reads = round(100 * n_reads / total_reads, 2),
        pct_umis = round(100 * n_umis / total_umis, 2)
    )]

    setorder(dt_major_freq, -n_reads)
    setnames(dt_major_freq, "major_status", "status")
    dt_major_freq[, roi := "major"]

    # === Minor ROI frequency table (if present) ===
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

            # Add percentages
            total_reads_minor <- sum(dt_minor_freq$n_reads)
            total_umis_minor <- sum(dt_minor_freq$n_umis)

            dt_minor_freq[, `:=`(
                pct_reads = round(100 * n_reads / total_reads_minor, 2),
                pct_umis = round(100 * n_umis / total_umis_minor, 2)
            )]

            setorder(dt_minor_freq, -n_reads)
            setnames(dt_minor_freq, "minor_status", "status")
            dt_minor_freq[, roi := "minor"]

            # Assess phasing viability
            # For phasing to work, we need both alleles (MUT and WT) present
            # with reasonable frequency (e.g., each >10% of UMIs)
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

    # Combine tables for easier viewing
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


# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Summarize SNV Detection Results
#'
#' Processes output from \code{\link{detect_snv}}, performs UMI error correction
#' using Levenshtein distance, and creates per-cell barcode summaries with read
#' counts. Does NOT perform mutation calling -- use \code{\link{flag_snv}} for
#' mutation calling with thresholds.
#'
#' @param session_name Character. Prefix for output files.
#' @param path_snv_table Character. Path to CSV output from \code{\link{detect_snv}}.
#' @param path_barcodes Character. Path to cell barcode whitelist (.csv or text file).
#' @param path_output_folder Character. Path to output directory for results.
#' @param umi_mismatch Integer. Maximum edit distance (Levenshtein) for UMI
#'   error correction. UMIs within this distance are collapsed to the most
#'   abundant sequence. Default: 2.
#' @param n_cores Integer. Number of CPU cores to use. Will be capped at
#'   \code{parallel::detectCores() - 1}. Default: 4.
#' @param skip_plots Logical. If \code{TRUE}, skip diagnostic plot generation.
#'   Default: \code{FALSE}.
#'
#' @return Invisibly returns a list with:
#' \describe{
#'   \item{summary}{data.table with per-CBC UMI counts (input for \code{flag_snv})}
#'   \item{umi_frequencies}{data.table with UMI read counts for thresholding}
#'   \item{knee_points}{List of detected knee/elbow/inflection points}
#'   \item{snp_frequencies}{List with major/minor ROI frequency tables}
#'   \item{phasing_viable}{Logical indicating if SNP phasing is viable}
#'   \item{phasing_message}{Character describing phasing assessment}
#'   \item{has_minor_roi}{Logical indicating if minor ROI data was found}
#'   \item{statistics}{List with QC statistics (n_input, n_valid, n_filtered, etc.)}
#' }
#'
#' @details
#' The output CSV uses semicolon (\code{;}) as separator with columns:
#' \code{cbc;umis;n_umi;n_umi_mut;n_umi_wt;n_snp;n_snp_mut;n_snp_wt}.
#' The \code{n_umi_mut} and \code{n_umi_wt} fields contain comma-separated
#' read counts per UMI.
#'
#' This output is designed as input for \code{\link{flag_snv}} which performs
#' mutation calling with customizable thresholds.
#'
#' @examples
#' \dontrun{
#' result <- summarize_snv(
#'   session_name = "EXP28",
#'   path_snv_table = "output/EXP28_snv_table.csv",
#'   path_barcodes = "barcodes/iPSC_final_barcodes.csv",
#'   path_output_folder = "output/",
#'   umi_mismatch = 2L,
#'   n_cores = 4L
#' )
#'
#' # View knee points for threshold selection
#' result$knee_points
#'
#' # Check phasing viability
#' message(result$phasing_message)
#' }
#'
#' @export
summarize_snv <- function(session_name,
                          path_snv_table,
                          path_barcodes,
                          path_output_folder,
                          umi_mismatch = 2L,
                          n_cores = 4L,
                          skip_plots = FALSE) {
    # =========================================================================
    # INPUT VALIDATION
    # =========================================================================

    stopifnot("`session_name` must be character" = is.character(session_name))
    stopifnot("`path_snv_table` not found" = file.exists(path_snv_table))
    stopifnot("`path_barcodes` not found" = file.exists(path_barcodes))
    stopifnot("`path_output_folder` not found" = dir.exists(path_output_folder))

    # Adjust cores
    available_cores <- parallel::detectCores()
    if (n_cores > available_cores - 1L) {
        n_cores <- max(1L, available_cores - 1L)
        message("n_cores adjusted to ", n_cores)
    }

    # =========================================================================
    # LOAD DATA
    # =========================================================================

    message(Sys.time(), " - Loading SNV table: ", path_snv_table)

    # Only load columns needed for summarization (not the full enriched table)
    all_cols <- names(fread(path_snv_table, nrows = 0))
    req_cols <- c("read_id", "CBC", "UMI", "major_base", "major_status")
    opt_cols <- c("minor_status")
    select_cols <- intersect(c(req_cols, opt_cols), all_cols)

    missing <- setdiff(req_cols, all_cols)
    if (length(missing) > 0) {
        stop("Missing required columns: ", paste(missing, collapse = ", "))
    }

    dt_input <- fread(path_snv_table, select = select_cols)
    n_input <- nrow(dt_input)
    message("  Loaded ", format(n_input, big.mark = ","), " reads (",
            length(select_cols), " of ", length(all_cols), " columns)")

    # Check for minor ROI (phasing data)
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

    # =========================================================================
    # FILTER BY BARCODES
    # =========================================================================

    message(Sys.time(), " - Loading barcode whitelist")

    # Handle both .csv and plain text formats
    if (grepl("\\.csv$", path_barcodes, ignore.case = TRUE)) {
        bc_raw <- fread(path_barcodes, header = FALSE)[[1]]
    } else {
        bc_raw <- readLines(path_barcodes)
    }

    # Extract 16bp barcode (handle -1 suffix from 10x)
    bc_16 <- substr(bc_raw, 1, 16)
    message("  Loaded ", length(bc_16), " barcodes")

    # Filter by valid major status and drop the full input table
    dt_valid <- dt_input[major_status %in% c("MUT", "WT")]
    n_valid <- nrow(dt_valid)
    rm(dt_input)
    gc(verbose = FALSE)
    message("  Valid major ROI calls: ", .fmt_pct(n_valid, n_input))

    # Filter by barcode whitelist
    setkey(dt_valid, CBC)
    dt_valid <- dt_valid[CBC %in% bc_16]
    n_filtered <- nrow(dt_valid)
    message("  After barcode filter: ", .fmt_pct(n_filtered, n_valid))
    message("  Unique CBCs: ", uniqueN(dt_valid$CBC))

    # =========================================================================
    # UMI CORRECTION
    # =========================================================================

    message(Sys.time(), " - Performing UMI correction (max mismatch = ", umi_mismatch, ")")

    # Keep essential columns only
    cols_keep <- c("read_id", "CBC", "UMI", "major_status")
    if (has_minor) cols_keep <- c(cols_keep, "minor_status")
    dt_valid <- dt_valid[, ..cols_keep]
    setnames(dt_valid, "UMI", "initial_UMI")

    # Update cols_keep to reflect renamed column
    cols_keep[cols_keep == "UMI"] <- "initial_UMI"

    n_cbcs <- uniqueN(dt_valid$CBC)
    dt_corrected <- dt_valid[,
        {
            n_initial_umis <- uniqueN(initial_UMI)
            result <- .correct_umis_optimized(.SD, umi_mismatch, n_cores)
            n_corrected_umis <- uniqueN(result$corrected_UMI)
            n_collapsed <- n_initial_umis - n_corrected_umis

            if (n_collapsed > 10) { # Log only significant corrections
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

    rm(dt_valid) # dt_input was already dropped after filtering
    gc(verbose = FALSE)

    # =========================================================================
    # COMPUTE UMI FREQUENCIES & KNEE POINTS
    # =========================================================================

    message(Sys.time(), " - Computing UMI frequencies (reads per UMI)")

    # Count reads per UMI (across all cells) - this is the duplication level
    umi_freq <- dt_corrected[, .(
        n_reads = .N
    ), by = .(CBC, corrected_UMI)]

    # Sort by reads per UMI (descending) and rank
    setorder(umi_freq, -n_reads)
    umi_freq[, rank := .I]

    message("  Total unique UMIs: ", nrow(umi_freq))
    message("  Max reads per UMI: ", max(umi_freq$n_reads))
    message("  Median reads per UMI: ", median(umi_freq$n_reads))

    # Detect knee points on reads per UMI distribution
    message(Sys.time(), " - Detecting knee/elbow/inflection points")
    knee_points <- .detect_knee_points(umi_freq$rank, umi_freq$n_reads)

    message("\n  === Knee Point Detection Results (Reads per UMI) ===")
    if (is.finite(knee_points$upper_knee)) message("  Upper Knee (start of steep drop): ", knee_points$upper_knee, " reads")
    if (is.finite(knee_points$inflection)) message("  Inflection Point (geometric): ", knee_points$inflection, " reads")
    if (is.finite(knee_points$lower_knee)) message("  Lower Knee (curve flattens): ", knee_points$lower_knee, " reads")
    if (is.finite(knee_points$uik)) message("  UIK (Unit Invariant Knee): ", knee_points$uik, " reads")
    if (is.finite(knee_points$cellranger)) message("  CellRanger-style (m/10): ", knee_points$cellranger, " reads")
    if (is.finite(knee_points$otsu)) message("  Otsu (inter-class variance): ", knee_points$otsu, " reads")
    if (is.finite(knee_points$mixture)) message("  Mixture model (EM intersection): ", knee_points$mixture, " reads")
    if (is.finite(knee_points$lower)) {
        message("  Suggested threshold range: ", knee_points$lower, " - ", knee_points$upper, " reads per UMI")
    }
    message("")

    # Note: Plots generated after dt_umi_counts is created (needed for 4th plot)

    # =========================================================================
    # SUMMARIZE PER CBC + UMI
    # =========================================================================

    message(Sys.time(), " - Summarizing reads per UMI")

    # Count MUT vs WT reads per corrected UMI
    dt_umi_counts <- dt_corrected[, .(
        n_mut = sum(major_status == "MUT"),
        n_wt = sum(major_status == "WT"),
        n_total = .N
    ), by = .(CBC, corrected_UMI)]

    # Add minor ROI counts if present
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

    # Generate diagnostic plots (moved here so dt_umi_counts is available for 4th plot)
    if (!skip_plots) {
        .plot_umi_distribution(umi_freq, knee_points, path_output_folder, session_name, dt_umi_counts)
    }

    # Sort by abundance within each CBC
    setorder(dt_umi_counts, CBC, -n_total)

    # Collapse to comma-separated strings per CBC
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

    # =========================================================================
    # SAVE OUTPUTS
    # =========================================================================

    message(Sys.time(), " - Saving outputs")

    # UMI summary (semicolon-separated - input for flag_snv)
    out_summary <- file.path(path_output_folder, paste0(session_name, "_snv_summary.csv"))
    fwrite(dt_summary, out_summary, sep = ";")
    message("  ", out_summary)

    # UMI frequencies (for thresholding decisions)
    out_freq <- file.path(path_output_folder, paste0(session_name, "_umi_frequencies.csv"))
    fwrite(umi_freq, out_freq)
    message("  ", out_freq)

    # Knee/threshold points
    knee_dt <- data.table(
        method = c("upper_knee", "inflection", "lower_knee", "uik",
                   "cellranger", "otsu", "mixture", "range_lower", "range_upper"),
        threshold = c(knee_points$upper_knee, knee_points$inflection,
                      knee_points$lower_knee, knee_points$uik,
                      knee_points$cellranger, knee_points$otsu,
                      knee_points$mixture, knee_points$lower, knee_points$upper)
    )
    out_knee <- file.path(path_output_folder, paste0(session_name, "_thresholds.csv"))
    fwrite(knee_dt, out_knee)
    message("  ", out_knee)

    # =========================================================================
    # CREATE AND SAVE SNP FREQUENCY TABLES
    # =========================================================================

    message(Sys.time(), " - Creating SNP frequency tables")

    snp_freq_tables <- .create_snp_frequency_tables(dt_corrected, has_minor)

    # Save major ROI frequency table
    out_major_freq <- file.path(path_output_folder, paste0(session_name, "_major_roi_frequency.csv"))
    fwrite(snp_freq_tables$major_freq, out_major_freq)
    message("  ", out_major_freq)

    # Save minor ROI frequency table if present
    if (!is.null(snp_freq_tables$minor_freq)) {
        out_minor_freq <- file.path(path_output_folder, paste0(session_name, "_minor_roi_frequency.csv"))
        fwrite(snp_freq_tables$minor_freq, out_minor_freq)
        message("  ", out_minor_freq)
    }

    # Save combined frequency table
    out_combined_freq <- file.path(path_output_folder, paste0(session_name, "_snp_frequency_combined.csv"))
    fwrite(snp_freq_tables$combined, out_combined_freq)
    message("  ", out_combined_freq)

    # Print phasing assessment
    message("\n  === SNP Frequency Summary ===")
    message("  Major ROI:")
    print(snp_freq_tables$major_freq[, .(status, n_reads, pct_reads, n_umis, pct_umis)])

    if (!is.null(snp_freq_tables$minor_freq)) {
        message("\n  Minor ROI:")
        print(snp_freq_tables$minor_freq[, .(status, n_reads, pct_reads, n_umis, pct_umis)])
    }

    message("\n  ", snp_freq_tables$phasing_message)

    # =========================================================================
    # RETURN RESULTS
    # =========================================================================

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
