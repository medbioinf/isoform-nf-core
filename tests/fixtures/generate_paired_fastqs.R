#!/usr/bin/env Rscript

# Rebuild the compact paired-end fixtures deterministically. Each 49 nt synthetic
# transcript yields an inward-facing 30 + 30 nt pair with a 19 nt overlap.

transcripts <- c(
    TX1 = "ATGCGTACGTTAGCTAGCTAACGATCGTACGATCGATCGTAGCTAGCTA",
    TX2 = "GCTAGCTAGGATCCGATCGTTACGATCGATCGTAGCTAACGTTAGCTAG",
    TX3 = "ACGATGCTACGTCGATGACCTAGCTAACGTGATCGTACGATCGTAGCTA",
    TX4 = "TGCATCGATGCTAGTCGATACGCTAGCTGATCGTACGATCGTAGCTAGC"
)

reverse_complement <- function(sequence) {
    paste(rev(chartr("ACGT", "TGCA", strsplit(sequence, "", fixed = TRUE)[[1]])), collapse = "")
}

sample_counts <- list(
    CONTROL_REP1 = c(TX1 = 240L, TX2 = 60L, TX3 = 150L, TX4 = 150L),
    CONTROL_REP2 = c(TX1 = 220L, TX2 = 80L, TX3 = 150L, TX4 = 150L),
    TREATED_REP1 = c(TX1 = 60L, TX2 = 240L, TX3 = 150L, TX4 = 150L),
    TREATED_REP2 = c(TX1 = 80L, TX2 = 220L, TX3 = 150L, TX4 = 150L)
)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
output_dir <- file.path(dirname(script_path), "reads")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_pairs <- function(sample_id, counts, r1_path, r2_path) {
    r1 <- gzfile(r1_path, open = "wt", compression = 9)
    r2 <- gzfile(r2_path, open = "wt", compression = 9)
    on.exit({
        close(r1)
        close(r2)
    }, add = TRUE)

    read_number <- 1L
    for (transcript_id in names(counts)) {
        sequence <- transcripts[[transcript_id]]
        read1 <- substr(sequence, 1L, 30L)
        read2 <- reverse_complement(substr(sequence, 20L, 49L))
        for (index in seq_len(counts[[transcript_id]])) {
            read_id <- sprintf("@%s_%s_%04d", tolower(sample_id), transcript_id, read_number)
            writeLines(c(read_id, read1, "+", strrep("F", 30L)), r1)
            writeLines(c(read_id, read2, "+", strrep("F", 30L)), r2)
            read_number <- read_number + 1L
        }
    }
}

for (sample_id in names(sample_counts)) {
    stem <- tolower(sample_id)
    write_pairs(
        sample_id,
        sample_counts[[sample_id]],
        file.path(output_dir, paste0(stem, "_paired_R1.fastq.gz")),
        file.path(output_dir, paste0(stem, "_paired_R2.fastq.gz"))
    )
}

# The lane fixtures exercise repeated-run grouping functions. Keep them valid
# gzip FASTQs as well as small enough for unit tests.
write_pairs(
    "CONTROL_REP1_L001",
    c(TX1 = 2L),
    file.path(output_dir, "control_rep1_L001_R1.fastq.gz"),
    file.path(output_dir, "control_rep1_L001_R2.fastq.gz")
)
write_pairs(
    "CONTROL_REP1_L002",
    c(TX2 = 2L),
    file.path(output_dir, "control_rep1_L002_R1.fastq.gz"),
    file.path(output_dir, "control_rep1_L002_R2.fastq.gz")
)
write_pairs(
    "TREATED_REP1_L001",
    c(TX2 = 2L),
    file.path(output_dir, "treated_rep1_L001_R1.fastq.gz"),
    file.path(output_dir, "treated_rep1_L001_R2.fastq.gz")
)
