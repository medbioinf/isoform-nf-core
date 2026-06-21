#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
    "Usage:",
    "Rscript run_deeploc2_prepare.R",
    "--isar-dir isar_analysis",
    "--outdir deeploc2_prepare",
    "[--top-n 25]",
    "[--qvalue-cutoff 0.05]",
    "[--dif-cutoff 0.1]",
    sep = "\n"
)

parse_args <- function(args) {
    opts <- list(
        top_n = "25",
        qvalue_cutoff = "0.05",
        dif_cutoff = "0.1"
    )
    i <- 1
    while (i <= length(args)) {
        key <- args[[i]]
        if (!startsWith(key, "--") || i == length(args)) {
            stop(usage, call. = FALSE)
        }
        opts[[sub("^--", "", key)]] <- args[[i + 1]]
        i <- i + 2
    }
    opts
}

first_existing <- function(paths) {
    hits <- paths[file.exists(paths)]
    if (length(hits) == 0) NULL else hits[[1]]
}

write_empty_csv <- function(path, columns) {
    empty <- as.data.frame(setNames(replicate(length(columns), character(0), simplify = FALSE), columns))
    write.csv(empty, path, row.names = FALSE)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

format_q <- function(q) {
    q <- safe_num(q)
    ifelse(is.na(q), "NA", formatC(q, format = "e", digits = 2))
}

format_dif_pp <- function(dif) {
    dif <- safe_num(dif)
    ifelse(is.na(dif), "NA", sprintf("%+.1f pp", dif * 100))
}

count_fasta_records <- function(path) {
    if (!file.exists(path)) {
        return(0L)
    }
    sum(grepl("^>", readLines(path, warn = FALSE)))
}

opts <- parse_args(args)
required <- c("isar-dir", "outdir")
missing <- required[!required %in% names(opts) | !nzchar(unlist(opts[required]))]
if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s\n%s", paste(missing, collapse = ", "), usage()), call. = FALSE)
}

isar_dir <- normalizePath(opts[["isar-dir"]], mustWork = TRUE)
outdir <- opts$outdir
top_n <- as.integer(opts$top_n)
qvalue_cutoff <- as.numeric(opts$qvalue_cutoff)
dif_cutoff <- as.numeric(opts$dif_cutoff)

if (is.na(top_n) || top_n < 1) {
    stop("--top-n must be a positive integer", call. = FALSE)
}
if (is.na(qvalue_cutoff) || qvalue_cutoff <= 0 || qvalue_cutoff >= 1) {
    stop("--qvalue-cutoff must be between 0 and 1", call. = FALSE)
}
if (is.na(dif_cutoff) || dif_cutoff < 0) {
    stop("--dif-cutoff must be a non-negative number", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rds_path <- first_existing(file.path(
    isar_dir,
    c("switchAnalyzeRlist_analyzed.rds", "switchAnalyzeRlist.rds", "switchAnalyzeRlist_imported.rds")
))

if (is.null(rds_path)) {
    write_empty_csv(
        file.path(outdir, "deeploc2_top_candidates.csv"),
        c("gene_id", "gene_name", "isoform_id", "condition_1", "condition_2", "IF1", "IF2", "dIF", "isoform_switch_q_value")
    )
    writeLines(
        c(
            "DeepLoc2 preparation summary",
            sprintf("Input ISAR directory: %s", isar_dir),
            "No switchAnalyzeRlist RDS file was found. No FASTA was generated."
        ),
        file.path(outdir, "deeploc2_prepare_notes.txt")
    )
    quit(save = "no", status = 0)
}

message("Loading switchAnalyzeRlist: ", rds_path)
switch_list <- readRDS(rds_path)
features <- switch_list$isoformFeatures

required_cols <- c("gene_id", "gene_name", "isoform_id", "condition_1", "condition_2", "IF1", "IF2", "dIF", "isoform_switch_q_value", "gene_switch_q_value")
missing_cols <- setdiff(required_cols, colnames(features))
if (length(missing_cols) > 0) {
    write_empty_csv(
        file.path(outdir, "deeploc2_top_candidates.csv"),
        c("gene_id", "gene_name", "isoform_id", "condition_1", "condition_2", "IF1", "IF2", "dIF", "isoform_switch_q_value")
    )
    writeLines(
        c(
            "DeepLoc2 preparation summary",
            sprintf("Input ISAR directory: %s", isar_dir),
            sprintf("Input RDS: %s", basename(rds_path)),
            sprintf("Skipped because isoformFeatures lacks switch-test columns: %s", paste(missing_cols, collapse = ", ")),
            "This is expected for import-only or no-switch runs."
        ),
        file.path(outdir, "deeploc2_prepare_notes.txt")
    )
    quit(save = "no", status = 0)
}

features$dIF <- safe_num(features$dIF)
features$isoform_switch_q_value <- safe_num(features$isoform_switch_q_value)
features$gene_switch_q_value <- safe_num(features$gene_switch_q_value)
features$is_deeploc2_candidate <- !is.na(features$isoform_switch_q_value) &
    features$isoform_switch_q_value < qvalue_cutoff &
    abs(features$dIF) >= dif_cutoff
features$direction <- ifelse(
    is.na(features$dIF) | features$dIF == 0,
    "unchanged",
    ifelse(features$dIF > 0, features$condition_2, features$condition_1)
)
features$dIF_pp <- format_dif_pp(features$dIF)
features$isoform_switch_q_value_formatted <- format_q(features$isoform_switch_q_value)
features$gene_switch_q_value_formatted <- format_q(features$gene_switch_q_value)

candidate_table <- features[features$is_deeploc2_candidate, , drop = FALSE]
candidate_table <- candidate_table[order(candidate_table$isoform_switch_q_value, -abs(candidate_table$dIF)), , drop = FALSE]
candidate_table <- utils::head(candidate_table, top_n)

candidate_cols <- intersect(
    c(
        "gene_id", "gene_name", "isoform_id", "condition_1", "condition_2",
        "IF1", "IF2", "dIF", "dIF_pp", "isoform_switch_q_value",
        "isoform_switch_q_value_formatted", "gene_switch_q_value",
        "gene_switch_q_value_formatted", "direction", "PTC", "iso_biotype"
    ),
    colnames(candidate_table)
)
utils::write.csv(
    candidate_table[, candidate_cols, drop = FALSE],
    file.path(outdir, "deeploc2_top_candidates.csv"),
    row.names = FALSE
)

output_prefix <- "isoform_deeploc2_candidates"
sequence_error <- tryCatch({
    switch_list_with_sequences <- extractSequence(
        switchAnalyzeRlist = switch_list,
        onlySwitchingGenes = TRUE,
        alpha = qvalue_cutoff,
        dIFcutoff = dif_cutoff,
        extractNTseq = FALSE,
        extractAAseq = TRUE,
        writeToFile = TRUE,
        pathToOutput = outdir,
        outputPrefix = output_prefix,
        quiet = TRUE
    )
    saveRDS(
        switch_list_with_sequences,
        file.path(outdir, "switchAnalyzeRlist_with_sequences.rds")
    )
    NULL
}, error = function(err) {
    conditionMessage(err)
})

aa_fasta <- file.path(outdir, paste0(output_prefix, "_AA.fasta"))
sequence_count <- count_fasta_records(aa_fasta)

notes <- c(
    "DeepLoc2 preparation summary",
    sprintf("Input ISAR directory: %s", isar_dir),
    sprintf("Input RDS: %s", basename(rds_path)),
    sprintf("q-value cutoff: %.3f", qvalue_cutoff),
    sprintf("dIF cutoff: %.3f", dif_cutoff),
    sprintf("Significant switch candidates in feature table: %d", sum(features$is_deeploc2_candidate, na.rm = TRUE)),
    sprintf("Top candidates written: %d", nrow(candidate_table)),
    sprintf("AA FASTA records written: %d", sequence_count),
    "",
    "Outputs:",
    "- isoform_deeploc2_candidates_AA.fasta: amino-acid FASTA for DeepLoc2.",
    "- deeploc2_top_candidates.csv: compact table of strongest isoforms.",
    "- switchAnalyzeRlist_with_sequences.rds: ISAR object after sequence extraction, when extraction succeeds."
)

if (!is.null(sequence_error)) {
    notes <- c(notes, "", "Sequence extraction warning:", sequence_error)
}

writeLines(notes, file.path(outdir, "deeploc2_prepare_notes.txt"))

message("Wrote DeepLoc2 preparation outputs to: ", outdir)
