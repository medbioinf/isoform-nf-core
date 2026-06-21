#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
    "Usage:",
    "Rscript run_pfam_import.R",
    "--isar-dir isar_analysis",
    "--pfam-results pfam_scan.out",
    "--outdir pfam_import",
    "[--qvalue-cutoff 0.05]",
    "[--dif-cutoff 0.1]",
    sep = "\n"
)

parse_args <- function(args) {
    opts <- list(qvalue_cutoff = "0.05", dif_cutoff = "0.1")
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

safe_num <- function(x) suppressWarnings(as.numeric(x))

format_q <- function(q) {
    q <- safe_num(q)
    ifelse(is.na(q), "NA", formatC(q, format = "e", digits = 2))
}

format_dif_pp <- function(dif) {
    dif <- safe_num(dif)
    ifelse(is.na(dif), "NA", sprintf("%+.1f pp", dif * 100))
}

opts <- parse_args(args)
required <- c("isar-dir", "pfam-results", "outdir")
missing <- required[!required %in% names(opts) | !nzchar(unlist(opts[required]))]
if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s\n%s", paste(missing, collapse = ", "), usage()), call. = FALSE)
}

isar_dir <- normalizePath(opts[["isar-dir"]], mustWork = TRUE)
pfam_result_file <- normalizePath(opts[["pfam-results"]], mustWork = TRUE)
outdir <- opts$outdir
qvalue_cutoff <- as.numeric(opts$qvalue_cutoff)
dif_cutoff <- as.numeric(opts$dif_cutoff)

if (is.na(qvalue_cutoff) || qvalue_cutoff <= 0 || qvalue_cutoff >= 1) {
    stop("--qvalue-cutoff must be between 0 and 1", call. = FALSE)
}
if (is.na(dif_cutoff) || dif_cutoff < 0) {
    stop("--dif-cutoff must be a non-negative number", call. = FALSE)
}

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

message("Importing Pfam results: ", pfam_result_file)
switch_list <- analyzePFAM(
    switchAnalyzeRlist = switch_list,
    pathToPFAMresultFile = pfam_result_file,
    showProgress = FALSE,
    quiet = TRUE
)

saveRDS(switch_list, file.path(outdir, "switchAnalyzeRlist_with_pfam.rds"))

if (!is.null(switch_list$domainAnalysis)) {
    utils::write.csv(
        switch_list$domainAnalysis,
        file.path(outdir, "pfam_domain_analysis.csv"),
        row.names = FALSE
    )
}

features <- switch_list$isoformFeatures
features$dIF <- safe_num(features$dIF)
features$isoform_switch_q_value <- safe_num(features$isoform_switch_q_value)
features$gene_switch_q_value <- safe_num(features$gene_switch_q_value)
features$is_significant_switch <- !is.na(features$isoform_switch_q_value) &
    features$isoform_switch_q_value < qvalue_cutoff &
    abs(features$dIF) >= dif_cutoff
features$dIF_pp <- format_dif_pp(features$dIF)
features$isoform_switch_q_value_formatted <- format_q(features$isoform_switch_q_value)
features$gene_switch_q_value_formatted <- format_q(features$gene_switch_q_value)

if (!"domain_identified" %in% colnames(features)) {
    features$domain_identified <- NA
}

feature_cols <- intersect(
    c(
        "gene_id", "gene_name", "isoform_id", "condition_1", "condition_2",
        "IF1", "IF2", "dIF", "dIF_pp", "isoform_switch_q_value",
        "isoform_switch_q_value_formatted", "gene_switch_q_value",
        "gene_switch_q_value_formatted", "domain_identified", "PTC", "iso_biotype"
    ),
    colnames(features)
)
utils::write.csv(
    features[, feature_cols, drop = FALSE],
    file.path(outdir, "isoform_features_with_pfam.csv"),
    row.names = FALSE
)

candidate_features <- features[features$is_significant_switch, , drop = FALSE]
candidate_features <- candidate_features[order(candidate_features$isoform_switch_q_value, -abs(candidate_features$dIF)), , drop = FALSE]
utils::write.csv(
    candidate_features[, feature_cols, drop = FALSE],
    file.path(outdir, "significant_switches_with_pfam.csv"),
    row.names = FALSE
)

domain_summary <- as.data.frame(table(
    ifelse(is.na(features$domain_identified), "unknown", as.character(features$domain_identified)),
    features$is_significant_switch
), stringsAsFactors = FALSE)
colnames(domain_summary) <- c("domain_identified", "significant_switch", "n_isoforms")
utils::write.csv(domain_summary, file.path(outdir, "pfam_domain_summary.csv"), row.names = FALSE)

domain_plot <- ggplot(domain_summary, aes(x = significant_switch, y = n_isoforms, fill = domain_identified)) +
    geom_col(position = "stack") +
    labs(
        title = "Pfam Domain Annotation Summary",
        subtitle = "Counts are grouped by significant switch status",
        x = "Significant isoform switch",
        y = "Number of isoforms",
        fill = "Domain identified"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

ggsave(file.path(outdir, "pfam_domain_summary.png"), domain_plot, width = 8, height = 5, dpi = 180)
ggsave(file.path(outdir, "pfam_domain_summary.pdf"), domain_plot, width = 8, height = 5)

consequence_error <- tryCatch({
    consequence_list <- analyzeSwitchConsequences(
        switchAnalyzeRlist = switch_list,
        consequencesToAnalyze = c("domains_identified"),
        alpha = qvalue_cutoff,
        dIFcutoff = dif_cutoff,
        removeNonConseqSwitches = FALSE,
        quiet = TRUE,
        showProgress = FALSE
    )
    if (!is.null(consequence_list$switchConsequence)) {
        utils::write.csv(
            consequence_list$switchConsequence,
            file.path(outdir, "pfam_switch_consequences.csv"),
            row.names = FALSE
        )
    }
    saveRDS(consequence_list, file.path(outdir, "switchAnalyzeRlist_with_pfam_consequences.rds"))
    NULL
}, error = function(err) {
    conditionMessage(err)
})

notes <- c(
    "Pfam import summary",
    sprintf("Input ISAR directory: %s", isar_dir),
    sprintf("Input RDS: %s", basename(rds_path)),
    sprintf("Pfam result file: %s", pfam_result_file),
    sprintf("q-value cutoff: %.3f", qvalue_cutoff),
    sprintf("dIF cutoff: %.3f", dif_cutoff),
    sprintf("Isoforms in feature table: %d", nrow(features)),
    sprintf("Significant switch candidates: %d", nrow(candidate_features)),
    sprintf("Isoforms with domain_identified == yes: %d", sum(features$domain_identified == "yes", na.rm = TRUE)),
    sprintf("Domain rows imported: %d", ifelse(is.null(switch_list$domainAnalysis), 0, nrow(switch_list$domainAnalysis)))
)

if (!is.null(consequence_error)) {
    notes <- c(notes, "Pfam consequence analysis warning:", consequence_error)
}

writeLines(notes, file.path(outdir, "pfam_import_notes.txt"))

message("Wrote Pfam import outputs to: ", outdir)
