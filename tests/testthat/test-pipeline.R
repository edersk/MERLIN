# End-to-end smoke test of the three-step pipeline on a tiny synthetic FASTQ that
# is generated on the fly, so no example data needs to ship with the package.
#
# Each read is laid out as the pipeline expects:
#   [adapter][cell barcode][UMI][poly(T)][major ROI][minor ROI][filler]
# with the ROI 'N' position filled with a chosen base (MUT or WT).

.make_read <- function(cbc, umi, major_base, minor_base) {
    adapter   <- "CTACACGACGCTCTTCCGATCT"
    poly_t    <- strrep("T", 12)
    filler    <- strrep("ACGT", 6)
    major_roi <- sub("N", major_base, "GAGATTTCNCTGTAGCT")
    minor_roi <- sub("N", minor_base, "AGAACAATNCCAAATGC")
    paste0(adapter, cbc, umi, poly_t, major_roi, minor_roi, filler)
}

test_that("detect_snv -> summarize_snv -> flag_snv runs end to end", {
    skip_on_cran()

    cbc_mut <- "AAAACCCCGGGGTTTT"
    cbc_wt  <- "ACACACACGTGTGTGT"

    # 12 identical reads per cell so each single UMI clears the read threshold.
    reads <- c(
        rep(.make_read(cbc_mut, "AAAAAAAAAAAA", "T", "C"), 12), # MUT cell
        rep(.make_read(cbc_wt,  "CCCCCCCCCCCC", "A", "T"), 12)  # WT cell
    )

    outdir <- file.path(tempdir(), paste0("merlin_test_", as.integer(Sys.time())))
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    on.exit(unlink(outdir, recursive = TRUE), add = TRUE)

    # Write the synthetic FASTQ and a barcode whitelist (with 10x -1 suffixes).
    fq <- file.path(outdir, "synth.fastq")
    qual <- vapply(reads, function(s) strrep("I", nchar(s)), character(1))
    writeLines(as.vector(rbind(paste0("@read", seq_along(reads)), reads, "+", qual)), fq)

    bc <- file.path(outdir, "barcodes.csv")
    writeLines(c(paste0(cbc_mut, "-1"), paste0(cbc_wt, "-1"), "GGGGAAAACCCCTTTT-1"), bc)

    # --- Step 1: detect ---
    detect_snv(
        session_name       = "test",
        path_input_file    = fq,
        path_output_folder = outdir,
        read1_seq          = "CTACACGACGCTCTTCCGATCT",
        read1_mismatch     = 2L,
        read1_within_range = 70L,
        major_roi_seq      = "GAGATTTCNCTGTAGCT",
        major_roi_mut      = "T", major_roi_wt = "A", major_roi_mismatch = 3L,
        minor_roi_seq      = "AGAACAATNCCAAATGC",
        minor_roi_mut      = "C", minor_roi_wt = "T", minor_roi_mismatch = 3L,
        len_CBC = 16L, len_UMI = 12L,
        n_reads_per_chunk = 1000L, n_cores = 1L
    )
    snv_table <- file.path(outdir, "test_snv_table.csv")
    expect_true(file.exists(snv_table))

    dt <- data.table::fread(snv_table)
    expect_equal(nrow(dt), 24L)
    expect_true(all(c("CBC", "UMI", "major_base", "major_status") %in% colnames(dt)))
    expect_true("MUT" %in% dt$major_status)
    expect_true("WT" %in% dt$major_status)

    # --- Step 2: summarize ---
    res <- summarize_snv(
        session_name       = "test",
        path_snv_table     = snv_table,
        path_barcodes      = bc,
        path_output_folder = outdir,
        umi_mismatch = 2L, n_cores = 1L, skip_plots = TRUE
    )
    expect_s3_class(res$summary, "data.table")
    expect_equal(nrow(res$summary), 2L)
    expect_true(all(c("cbc", "n_umi_mut", "n_umi_wt") %in% colnames(res$summary)))

    # --- Step 3: flag ---
    flags <- flag_snv(
        path_snv_summary   = res$summary,
        session_name       = "test",
        threshold = 10, consensus_threshold = 0.8,
        use_phasing = TRUE, minor_percentage = 0.2, minor_on_major = "MUT"
    )
    expect_s3_class(flags, "data.table")
    expect_equal(nrow(flags), 2L)
    expect_true(all(c(
        "cbc", "barcode_10x", "mutation_call", "phased_call", "phased_call_simple"
    ) %in% colnames(flags)))

    # The MUT cell is called MUT (and phased MUT); the WT cell is called WT.
    call_mut <- flags$mutation_call[flags$cbc == cbc_mut]
    call_wt  <- flags$mutation_call[flags$cbc == cbc_wt]
    expect_equal(call_mut, "MUT")
    expect_equal(call_wt, "WT")
    expect_equal(flags$phased_call_simple[flags$cbc == cbc_mut], "MUT")
})
