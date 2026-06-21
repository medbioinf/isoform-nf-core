#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
    "Usage:",
    "Rscript run_iupred2a_import.R",
    "--isar-dir isar_analysis",
    "--iupred2a-results iupred2a_anchor2_isar.out",
    "--outdir iupred2a_import",
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
required <- c("isar-dir", "iupred2a-results", "outdir")
missing <- required[!required %in% names(opts) | !nzchar(unlist(opts[required]))]
if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s\n%s", paste(missing, collapse = ", "), usage()), call. = FALSE)
}

isar_dir <- normalizePath(opts[["isar-dir"]], mustWork = TRUE)
iupred2a_result_file <- normalizePath(opts[["iupred2a-results"]], mustWork = TRUE)
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

message("Importing IUPred2A/ANCHOR2 results: ", iupred2a_result_file)
switch_list <- analyzeIUPred2A(
    switchAnalyzeRlist = switch_list,
    pathToIUPred2AresultFile = iupred2a_result_file,
    showProgress = FALSE
)

saveRDS(switch_list, file.path(outdir, "switchAnalyzeRlist_with_iupred2a.rds"))

idr_rows <- if (!is.null(switch_list$idrAnalysis)) switch_list$idrAnalysis else data.frame()
if (nrow(idr_rows) > 0) {
    utils::write.csv(
        idr_rows,
        file.path(outdir, "iupred2a_idr_analysis.csv"),
        row.names = FALSE
    )
} else {
    write_empty_csv(file.path(outdir, "iupred2a_idr_analysis.csv"), c("isoform_id"))
}

features <- switch_list$isoformFeatures
feature_cols <- intersect(
    c(
        "gene_id", "gene_name", "isoform_id", "condition_1", "condition_2",
        "IF1", "IF2", "dIF", "isoform_switch_q_value", "gene_switch_q_value",
        "IDR_identified", "IDR_type", "PTC", "iso_biotype"
    ),
    colnames(features)
)
utils::write.csv(
    features[, feature_cols, drop = FALSE],
    file.path(outdir, "isoform_features_with_iupred2a.csv"),
    row.names = FALSE
)

if ("IDR_identified" %in% colnames(features)) {
    idr_summary <- as.data.frame(table(
        ifelse(is.na(features$IDR_identified), "unknown", as.character(features$IDR_identified))
    ), stringsAsFactors = FALSE)
    colnames(idr_summary) <- c("IDR_identified", "n_isoforms")
} else {
    idr_summary <- data.frame(IDR_identified = "not_available", n_isoforms = nrow(features))
}
utils::write.csv(idr_summary, file.path(outdir, "iupred2a_idr_summary.csv"), row.names = FALSE)

idr_plot <- ggplot(idr_summary, aes(x = IDR_identified, y = n_isoforms, fill = IDR_identified)) +
    geom_col(width = 0.72) +
    labs(
        title = "IUPred2A Disorder Annotation Summary",
        subtitle = "Isoforms grouped by imported IDR annotation availability",
        x = "IDR identified",
        y = "Number of isoforms",
        fill = "IDR identified"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

ggsave(file.path(outdir, "iupred2a_idr_summary.png"), idr_plot, width = 7, height = 4.5, dpi = 180)
ggsave(file.path(outdir, "iupred2a_idr_summary.pdf"), idr_plot, width = 7, height = 4.5)

notes <- c(
    "IUPred2A import summary",
    sprintf("Input ISAR directory: %s", isar_dir),
    sprintf("Input RDS: %s", basename(rds_path)),
    sprintf("IUPred2A result file: %s", iupred2a_result_file),
    sprintf("Isoforms in feature table: %d", nrow(features)),
    sprintf("IDR rows imported: %d", nrow(idr_rows)),
    sprintf("Isoforms with IDR_identified == yes: %d", ifelse("IDR_identified" %in% colnames(features), sum(features$IDR_identified == "yes", na.rm = TRUE), 0)),
    "",
    "Interpretation:",
    "IUPred2A predicts intrinsically disordered protein regions.",
    "ANCHOR2 adds predicted disordered binding regions.",
    "These annotations can be rendered by IsoformSwitchAnalyzeR::switchPlot() as IDR tracks."
)

writeLines(notes, file.path(outdir, "iupred2a_import_notes.txt"))

message("Wrote IUPred2A import outputs to: ", outdir)
