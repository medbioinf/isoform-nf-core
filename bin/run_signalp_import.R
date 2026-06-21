#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
    "Usage:",
    "Rscript run_signalp_import.R",
    "--isar-dir isar_analysis",
    "--signalp-results signalp5_summary.signalp5",
    "--outdir signalp_import",
    sep = "\n"
)

parse_args <- function(args) {
    opts <- list()
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

opts <- parse_args(args)
required <- c("isar-dir", "signalp-results", "outdir")
missing <- required[!required %in% names(opts) | !nzchar(unlist(opts[required]))]
if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s\n%s", paste(missing, collapse = ", "), usage()), call. = FALSE)
}

isar_dir <- normalizePath(opts[["isar-dir"]], mustWork = TRUE)
signalp_result_file <- normalizePath(opts[["signalp-results"]], mustWork = TRUE)
outdir <- opts$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rds_path <- first_existing(file.path(
    isar_dir,
    c("switchAnalyzeRlist_analyzed.rds", "switchAnalyzeRlist.rds", "switchAnalyzeRlist_imported.rds")
))

if (is.null(rds_path)) {
    stop("No switchAnalyzeRlist RDS file found in the ISAR result directory.", call. = FALSE)
}

message("Loading switchAnalyzeRlist: ", rds_path)
switch_list <- readRDS(rds_path)

message("Importing SignalP results: ", signalp_result_file)
switch_list <- analyzeSignalP(
    switchAnalyzeRlist = switch_list,
    pathToSignalPresultFile = signalp_result_file
)

saveRDS(switch_list, file.path(outdir, "switchAnalyzeRlist_with_signalp.rds"))

signalp_rows <- if (!is.null(switch_list$signalPeptideAnalysis)) switch_list$signalPeptideAnalysis else data.frame()
if (nrow(signalp_rows) > 0) {
    utils::write.csv(
        signalp_rows,
        file.path(outdir, "signalp_signal_peptide_analysis.csv"),
        row.names = FALSE
    )
} else {
    write_empty_csv(file.path(outdir, "signalp_signal_peptide_analysis.csv"), c("isoform_id"))
}

features <- switch_list$isoformFeatures
feature_cols <- intersect(
    c(
        "gene_id", "gene_name", "isoform_id", "condition_1", "condition_2",
        "IF1", "IF2", "dIF", "isoform_switch_q_value", "gene_switch_q_value",
        "signal_peptide_identified", "PTC", "iso_biotype"
    ),
    colnames(features)
)
utils::write.csv(
    features[, feature_cols, drop = FALSE],
    file.path(outdir, "isoform_features_with_signalp.csv"),
    row.names = FALSE
)

if ("signal_peptide_identified" %in% colnames(features)) {
    signalp_summary <- as.data.frame(table(
        ifelse(is.na(features$signal_peptide_identified), "unknown", as.character(features$signal_peptide_identified))
    ), stringsAsFactors = FALSE)
    colnames(signalp_summary) <- c("signal_peptide_identified", "n_isoforms")
} else {
    signalp_summary <- data.frame(signal_peptide_identified = "not_available", n_isoforms = nrow(features))
}
utils::write.csv(signalp_summary, file.path(outdir, "signalp_summary.csv"), row.names = FALSE)

signalp_plot <- ggplot(signalp_summary, aes(x = signal_peptide_identified, y = n_isoforms, fill = signal_peptide_identified)) +
    geom_col(width = 0.72) +
    labs(
        title = "SignalP Signal Peptide Annotation Summary",
        subtitle = "Isoforms grouped by imported signal peptide annotation availability",
        x = "Signal peptide identified",
        y = "Number of isoforms",
        fill = "Signal peptide identified"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

ggsave(file.path(outdir, "signalp_summary.png"), signalp_plot, width = 7, height = 4.5, dpi = 180)
ggsave(file.path(outdir, "signalp_summary.pdf"), signalp_plot, width = 7, height = 4.5)

notes <- c(
    "SignalP import summary",
    sprintf("Input ISAR directory: %s", isar_dir),
    sprintf("Input RDS: %s", basename(rds_path)),
    sprintf("SignalP result file: %s", signalp_result_file),
    sprintf("Isoforms in feature table: %d", nrow(features)),
    sprintf("Signal peptide rows imported: %d", nrow(signalp_rows)),
    sprintf("Isoforms with signal_peptide_identified == yes: %d", ifelse("signal_peptide_identified" %in% colnames(features), sum(features$signal_peptide_identified == "yes", na.rm = TRUE), 0)),
    "",
    "Interpretation:",
    "SignalP predicts signal peptides: short protein segments that can route proteins into the secretory pathway.",
    "These annotations can be rendered by IsoformSwitchAnalyzeR::switchPlot() as signal peptide tracks."
)

writeLines(notes, file.path(outdir, "signalp_import_notes.txt"))

message("Wrote SignalP import outputs to: ", outdir)
