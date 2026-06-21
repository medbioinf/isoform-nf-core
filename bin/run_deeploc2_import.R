#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
    "Usage:",
    "Rscript run_deeploc2_import.R",
    "--isar-dir isar_analysis",
    "--deeploc2-results deeploc2_isar.csv",
    "--outdir deeploc2_import",
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
    utils::write.csv(empty, path, row.names = FALSE)
}

opts <- parse_args(args)
required <- c("isar-dir", "deeploc2-results", "outdir")
missing <- required[!required %in% names(opts) | !nzchar(unlist(opts[required]))]
if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s\n%s", paste(missing, collapse = ", "), usage()), call. = FALSE)
}

isar_dir <- normalizePath(opts[["isar-dir"]], mustWork = TRUE)
deeploc2_result_file <- normalizePath(opts[["deeploc2-results"]], mustWork = TRUE)
outdir <- opts$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rds_path <- first_existing(file.path(
    isar_dir,
    c(
        "switchAnalyzeRlist_with_deeploc2.rds",
        "switchAnalyzeRlist_with_deeptmhmm.rds",
        "switchAnalyzeRlist_with_signalp.rds",
        "switchAnalyzeRlist_with_iupred2a.rds",
        "switchAnalyzeRlist_with_pfam_consequences.rds",
        "switchAnalyzeRlist_with_pfam.rds",
        "switchAnalyzeRlist_analyzed.rds",
        "switchAnalyzeRlist.rds",
        "switchAnalyzeRlist_imported.rds"
    )
))

if (is.null(rds_path)) {
    stop("No switchAnalyzeRlist RDS file found in the ISAR result directory.", call. = FALSE)
}

message("Loading switchAnalyzeRlist: ", rds_path)
switch_list <- readRDS(rds_path)

message("Importing DeepLoc2 results: ", deeploc2_result_file)
switch_list <- analyzeDeepLoc2(
    switchAnalyzeRlist = switch_list,
    pathToDeepLoc2resultFile = deeploc2_result_file
)

saveRDS(switch_list, file.path(outdir, "switchAnalyzeRlist_with_deeploc2.rds"))

features <- switch_list$isoformFeatures
feature_cols <- intersect(
    c(
        "gene_id", "gene_name", "isoform_id", "condition_1", "condition_2",
        "IF1", "IF2", "dIF", "isoform_switch_q_value", "gene_switch_q_value",
        "sub_cell_location", "PTC", "iso_biotype"
    ),
    colnames(features)
)
utils::write.csv(
    features[, feature_cols, drop = FALSE],
    file.path(outdir, "isoform_features_with_deeploc2.csv"),
    row.names = FALSE
)

location_rows <- if ("sub_cell_location" %in% colnames(features)) {
    features[
        !is.na(features$sub_cell_location) & features$sub_cell_location != "",
        intersect(c("gene_id", "gene_name", "isoform_id", "sub_cell_location"), colnames(features)),
        drop = FALSE
    ]
} else {
    data.frame()
}
if (nrow(location_rows) > 0) {
    utils::write.csv(location_rows, file.path(outdir, "deeploc2_location_analysis.csv"), row.names = FALSE)
} else {
    write_empty_csv(file.path(outdir, "deeploc2_location_analysis.csv"), c("gene_id", "gene_name", "isoform_id", "sub_cell_location"))
}

if ("sub_cell_location" %in% colnames(features)) {
    location_values <- ifelse(is.na(features$sub_cell_location) | features$sub_cell_location == "", "unknown", as.character(features$sub_cell_location))
    location_summary <- as.data.frame(table(location_values), stringsAsFactors = FALSE)
    colnames(location_summary) <- c("sub_cell_location", "n_isoforms")
} else {
    location_summary <- data.frame(sub_cell_location = "not_available", n_isoforms = nrow(features))
}
location_summary <- location_summary[order(-location_summary$n_isoforms, location_summary$sub_cell_location), , drop = FALSE]
utils::write.csv(location_summary, file.path(outdir, "deeploc2_summary.csv"), row.names = FALSE)

summary_plot <- ggplot(utils::head(location_summary, 12), aes(x = reorder(sub_cell_location, n_isoforms), y = n_isoforms, fill = sub_cell_location)) +
    geom_col(width = 0.72) +
    coord_flip() +
    labs(
        title = "DeepLoc2 Subcellular Location Summary",
        subtitle = "Top imported location labels by isoform count",
        x = "Subcellular location",
        y = "Number of isoforms",
        fill = "Location"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

ggsave(file.path(outdir, "deeploc2_summary.png"), summary_plot, width = 7, height = 4.5, dpi = 180)
ggsave(file.path(outdir, "deeploc2_summary.pdf"), summary_plot, width = 7, height = 4.5)

notes <- c(
    "DeepLoc2 import summary",
    sprintf("Input ISAR directory: %s", isar_dir),
    sprintf("Input RDS: %s", basename(rds_path)),
    sprintf("DeepLoc2 result file: %s", deeploc2_result_file),
    sprintf("Isoforms in feature table: %d", nrow(features)),
    sprintf("Isoforms with sub_cell_location: %d", ifelse("sub_cell_location" %in% colnames(features), sum(!is.na(features$sub_cell_location) & features$sub_cell_location != ""), 0)),
    "",
    "Interpretation:",
    "DeepLoc2 predicts subcellular localization. These annotations can be rendered by IsoformSwitchAnalyzeR::switchPlot() as localization tracks."
)

writeLines(notes, file.path(outdir, "deeploc2_import_notes.txt"))

message("Wrote DeepLoc2 import outputs to: ", outdir)
