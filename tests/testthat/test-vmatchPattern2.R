test_that("vmatchPattern2 returns one IRanges per subject sequence", {
    x <- Biostrings::DNAStringSet(c("AAGCGCGATATG", "GCAAAATCCCC"))
    hits <- vmatchPattern2(Biostrings::DNAString("GATATG"), x, max.mismatch = 0)
    expect_length(hits, 2L)
})

test_that("vmatchPattern2 finds exact matches", {
    x <- Biostrings::DNAStringSet(c("AAGCGCGATATG", "GCAAAATCCCC"))

    # "GATATG" occurs once in the first sequence and not at all in the second.
    hits <- vmatchPattern2(Biostrings::DNAString("GATATG"), x, max.mismatch = 0)
    expect_equal(length(hits[[1]]), 1L)
    expect_equal(length(hits[[2]]), 0L)

    # Overlapping matches are all reported: "GCG" hits twice in "AAGCGCGATATG".
    overlap <- vmatchPattern2(Biostrings::DNAString("GCG"), x, max.mismatch = 0)
    expect_equal(length(overlap[[1]]), 2L)
})

test_that("vmatchPattern2 honours the mismatch tolerance", {
    x <- Biostrings::DNAStringSet(c("AAGCGCGATATG", "GCAAAATCCCC"))

    # "GATTTG" differs from "GATATG" by a single base.
    none <- vmatchPattern2(Biostrings::DNAString("GATTTG"), x, max.mismatch = 0)
    expect_equal(length(none[[1]]), 0L)

    one <- vmatchPattern2(Biostrings::DNAString("GATTTG"), x, max.mismatch = 1)
    expect_equal(length(one[[1]]), 1L)
})
