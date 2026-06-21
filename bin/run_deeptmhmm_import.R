#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
    "Usage:",
    "Rscript run_deeptmhmm_import.R",
    "--isar-dir isar_analysis",
    "--deeptmhmm-results deeptmhmm_regions_isar.tsv",
    "--outdir deeptmhmm_import",
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
required <- c("isar-dir", "deeptmhmm-results", "outdir")
missing <- required[!required %in% names(opts) | !nzchar(unlist(opts[required]))]
if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s\n%s", paste(missing, collapse = ", "), usage()), call. = FALSE)
}

isar_dir <- normalizePath(opts[["isar-dir"]], mustWork = TRUE)
deeptmhmm_result_file <- normalizePath(opts[["deeptmhmm-results"]], mustWork = TRUE)
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

message("Importing DeepTMHMM results: ", deeptmhmm_result_file)
switch_list <- analyzeDeepTMHMM(
    switchAnalyzeRlist = switch_list,
    pathToDeepTMHMMresultFile = deeptmhmm_result_file
)

saveRDS(switch_list, file.path(outdir, "switchAnalyzeRlist_with_deeptmhmm.rds"))

topology_rows <- if (!is.null(switch_list$topologyAnalysis)) switch_list$topologyAnalysis else data.frame()
if (nrow(topology_rows) > 0) {
    utils::write.csv(topology_rows, file.path(outdir, "deeptmhmm_topology_analysis.csv"), row.names = FALSE)
} else {
    write_empty_csv(file.path(outdir, "deeptmhmm_topology_analysis.csv"), c("isoform_id", "topology", "start", "end"))
}

features <- switch_list$isoformFeatures
feature_cols <- intersect(
    c(
        "gene_id", "gene_name", "isoform_id", "condition_1", "condition_2",
        "IF1", "IF2", "dIF", "isoform_switch_q_value", "gene_switch_q_value",
        "topology_identified", "PTC", "iso_biotype"
    ),
    colnames(features)
)
utils::write.csv(
    features[, feature_cols, drop = FALSE],
    file.path(outdir, "isoform_features_with_deeptmhmm.csv"),
    row.names = FALSE
)

if ("topology_identified" %in% colnames(features)) {
    topology_summary <- as.data.frame(table(
        ifelse(is.na(features$topology_identified), "unknown", as.character(features$topology_identified))
    ), stringsAsFactors = FALSE)
    colnames(topology_summary) <- c("topology_identified", "n_isoforms")
} else {
    topology_summary <- data.frame(topology_identified = "not_available", n_isoforms = nrow(features))
}
utils::write.csv(topology_summary, file.path(outdir, "deeptmhmm_summary.csv"), row.names = FALSE)

summary_plot <- ggplot(topology_summary, aes(x = topology_identified, y = n_isoforms, fill = topology_identified)) +
    geom_col(width = 0.72) +
    labs(
        title = "DeepTMHMM Topology Annotation Summary",
        subtitle = "Isoforms grouped by imported topology annotation availability",
        x = "Topology identified",
        y = "Number of isoforms",
        fill = "Topology identified"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")

ggsave(file.path(outdir, "deeptmhmm_summary.png"), summary_plot, width = 7, height = 4.5, dpi = 180)
ggsave(file.path(outdir, "deeptmhmm_summary.pdf"), summary_plot, width = 7, height = 4.5)

notes <- c(
    "DeepTMHMM import summary",
    sprintf("Input ISAR directory: %s", isar_dir),
    sprintf("Input RDS: %s", basename(rds_path)),
    sprintf("DeepTMHMM result file: %s", deeptmhmm_result_file),
    sprintf("Isoforms in feature table: %d", nrow(features)),
    sprintf("Topology rows imported: %d", nrow(topology_rows)),
    sprintf("Isoforms with topology_identified == yes: %d", ifelse("topology_identified" %in% colnames(features), sum(features$topology_identified == "yes", na.rm = TRUE), 0)),
    "",
    "Interpretation:",
    "DeepTMHMM predicts transmembrane topology regions. These annotations can be rendered by IsoformSwitchAnalyzeR::switchPlot() as topology tracks."
)

writeLines(notes, file.path(outdir, "deeptmhmm_import_notes.txt"))

message("Wrote DeepTMHMM import outputs to: ", outdir)
