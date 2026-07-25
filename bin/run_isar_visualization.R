#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
    "Usage:",
    "Rscript run_isar_visualization.R",
    "--isar-dir isar_analysis",
    "--outdir isar_visualization",
    "[--top-n 10]",
    "[--qvalue-cutoff 0.05]",
    "[--dif-cutoff 0.1]",
    "[--condition-1 control --condition-2 case --comparison-label case_vs_control]",
    sep = "\n"
)

parse_args <- function(args) {
    opts <- list(
        top_n = "10",
        qvalue_cutoff = "0.05",
        dif_cutoff = "0.1",
        condition_1 = "",
        condition_2 = "",
        comparison_label = ""
    )
    i <- 1
    while (i <= length(args)) {
        key <- args[[i]]
        if (!startsWith(key, "--") || i == length(args)) {
            stop(usage, call. = FALSE)
        }
        option_name <- gsub("-", "_", sub("^--", "", key), fixed = TRUE)
        opts[[option_name]] <- args[[i + 1]]
        i <- i + 2
    }
    opts
}

write_empty_csv <- function(path, columns) {
    empty <- as.data.frame(setNames(replicate(length(columns), character(0), simplify = FALSE), columns))
    write.csv(empty, path, row.names = FALSE)
}

write_text <- function(path, lines) {
    writeLines(as.character(lines), path)
}

first_existing <- function(paths) {
    hits <- paths[file.exists(paths)]
    if (length(hits) == 0) NULL else hits[[1]]
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

safe_neg_log10 <- function(x) {
    x <- safe_num(x)
    finite_positive <- x[is.finite(x) & x > 0]
    floor_value <- if (length(finite_positive) > 0) min(finite_positive) / 10 else .Machine$double.xmin
    -log10(pmax(x, floor_value, na.rm = TRUE))
}

q_label <- function(q, qvalue_cutoff) {
    q <- safe_num(q)
    ifelse(is.na(q), "NA", ifelse(q <= qvalue_cutoff, "significant", "ns"))
}

format_q <- function(q) {
    q <- safe_num(q)
    ifelse(is.na(q), "NA", formatC(q, format = "e", digits = 2))
}

format_dif_pp <- function(dif) {
    dif <- safe_num(dif)
    ifelse(is.na(dif), "NA", sprintf("%+.1f pp", dif * 100))
}

format_if_pct <- function(value) {
    value <- safe_num(value)
    ifelse(is.na(value), "NA", sprintf("%.1f%%", value * 100))
}

short_isoform_id <- function(ids) sub("\\.[0-9]+$", "", as.character(ids))

sanitize_filename <- function(value) {
    result <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(value))
    result[result == ""] <- "comparison"
    result
}

make_switch_key <- function(gene_id, condition_1, condition_2) {
    paste(as.character(gene_id), as.character(condition_1), as.character(condition_2), sep = "||")
}

make_comparison_label <- function(condition_1, condition_2) {
    sprintf("%s vs %s", as.character(condition_2), as.character(condition_1))
}

script_path <- function() {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg) != 1) {
        stop("Could not determine the visualization script path.", call. = FALSE)
    }
    normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}

comparison_table <- function(features, isar_dir) {
    pairs <- unique(features[, c("condition_1", "condition_2"), drop = FALSE])
    pairs <- pairs[stats::complete.cases(pairs) & nzchar(pairs$condition_1) & nzchar(pairs$condition_2), , drop = FALSE]
    pairs$comparison <- sprintf("%s_vs_%s", pairs$condition_2, pairs$condition_1)

    comparison_file <- file.path(isar_dir, "comparisons.csv")
    if (file.exists(comparison_file)) {
        labels <- read.csv(comparison_file, stringsAsFactors = FALSE, check.names = FALSE)
        required <- c("contrast", "condition_1", "condition_2")
        if (all(required %in% colnames(labels))) {
            for (i in seq_len(nrow(pairs))) {
                match_row <- labels[
                    labels$condition_1 == pairs$condition_1[[i]] &
                        labels$condition_2 == pairs$condition_2[[i]],
                    ,
                    drop = FALSE
                ]
                if (nrow(match_row) > 0 && nzchar(match_row$contrast[[1]])) {
                    pairs$comparison[[i]] <- match_row$contrast[[1]]
                }
            }
        }
    }

    pairs$output_directory <- make.unique(sanitize_filename(pairs$comparison), sep = "_")
    pairs
}

dispatch_comparisons <- function(comparisons, opts, isar_dir, outdir) {
    statuses <- integer(nrow(comparisons))
    rscript <- file.path(R.home("bin"), "Rscript")
    own_script <- script_path()

    for (i in seq_len(nrow(comparisons))) {
        comparison_outdir <- file.path(outdir, comparisons$output_directory[[i]])
        child_args <- c(
            own_script,
            "--isar-dir", isar_dir,
            "--outdir", comparison_outdir,
            "--top-n", opts$top_n,
            "--qvalue-cutoff", opts$qvalue_cutoff,
            "--dif-cutoff", opts$dif_cutoff,
            "--condition-1", comparisons$condition_1[[i]],
            "--condition-2", comparisons$condition_2[[i]],
            "--comparison-label", comparisons$comparison[[i]]
        )
        statuses[[i]] <- system2(rscript, args = shQuote(child_args))
    }

    comparisons$status <- ifelse(statuses == 0, "complete", "failed")
    write.csv(comparisons, file.path(outdir, "comparison_visualizations.csv"), row.names = FALSE)
    writeLines(
        c(
            "ISAR per-comparison visualization summary",
            sprintf("Comparisons discovered: %d", nrow(comparisons)),
            sprintf("Comparisons completed: %d", sum(statuses == 0)),
            "Each comparison directory contains an independently ranked top-N set and its own volcano plot."
        ),
        file.path(outdir, "visualization_notes.txt")
    )

    if (any(statuses != 0)) {
        stop("One or more per-comparison visualizations failed. See comparison_visualizations.csv.", call. = FALSE)
    }
}

write_no_switch_outputs <- function(outdir, notes) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    write_empty_csv(
        file.path(outdir, "top_isoform_candidates.csv"),
        c("gene_id", "gene_name", "isoform_id", "condition_1", "condition_2", "IF1", "IF2", "dIF", "isoform_switch_q_value", "gene_switch_q_value")
    )
    write_empty_csv(
        file.path(outdir, "top_gene_summary.csv"),
        c("gene_id", "gene_name", "condition_1", "condition_2", "gene_switch_q_value", "isoforms_in_gene", "significant_switch_isoforms", "max_abs_dIF")
    )
    write_empty_csv(
        file.path(outdir, "ptc_switch_summary.csv"),
        c("PTC_status", "direction", "n")
    )
    writeLines(notes, file.path(outdir, "visualization_notes.txt"))
}

opts <- parse_args(args)
required <- c("isar_dir", "outdir")
missing <- required[!required %in% names(opts) | !nzchar(unlist(opts[required]))]
if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s\n%s", paste(missing, collapse = ", "), usage()), call. = FALSE)
}

isar_dir <- normalizePath(opts$isar_dir, mustWork = TRUE)
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
if (is.na(dif_cutoff) || dif_cutoff < 0 || dif_cutoff > 1) {
    stop("--dif-cutoff must be between 0 and 1 inclusive", call. = FALSE)
}
has_condition_1 <- nzchar(opts$condition_1)
has_condition_2 <- nzchar(opts$condition_2)
if (xor(has_condition_1, has_condition_2)) {
    stop("--condition-1 and --condition-2 must be supplied together.", call. = FALSE)
}
cutoff_description <- sprintf("q <= %s and |dIF| >= %s", format(qvalue_cutoff), format(dif_cutoff))

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

notes <- c(
    "ISAR visualization summary",
    sprintf("Input directory: %s", isar_dir)
)

analyzed_rds <- file.path(isar_dir, "switchAnalyzeRlist_analyzed.rds")
imported_rds <- file.path(isar_dir, "switchAnalyzeRlist_imported.rds")
rds_path <- first_existing(c(analyzed_rds, imported_rds))

if (is.null(rds_path)) {
    write_no_switch_outputs(outdir, c(notes, "No switchAnalyzeRlist RDS file was found. No plots were generated."))
    quit(save = "no", status = 0)
}

if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)) {
    stop("IsoformSwitchAnalyzeR is required to read the switchAnalyzeRlist RDS.", call. = FALSE)
}

switch_list <- readRDS(rds_path)
features <- switch_list$isoformFeatures

if (is.null(features) || !is.data.frame(features) || nrow(features) == 0) {
    write_no_switch_outputs(outdir, c(notes, sprintf("RDS file: %s", basename(rds_path)), "No isoformFeatures table was found. No plots were generated."))
    quit(save = "no", status = 0)
}

required_feature_cols <- c("gene_id", "gene_name", "isoform_id", "condition_1", "condition_2", "IF1", "IF2", "dIF", "isoform_switch_q_value", "gene_switch_q_value")
missing_feature_cols <- setdiff(required_feature_cols, colnames(features))
if (length(missing_feature_cols) > 0) {
    write_no_switch_outputs(
        outdir,
        c(
            notes,
            sprintf("RDS file: %s", basename(rds_path)),
            sprintf("Visualization skipped because isoformFeatures lacks switch-test columns: %s", paste(missing_feature_cols, collapse = ", ")),
            "This is expected for import-only or no-switch runs."
        )
    )
    quit(save = "no", status = 0)
}

if (!has_condition_1) {
    comparisons <- comparison_table(features, isar_dir)
    if (nrow(comparisons) == 0) {
        write_no_switch_outputs(
            outdir,
            c(notes, sprintf("RDS file: %s", basename(rds_path)), "No tested comparisons were found. No plots were generated.")
        )
        quit(save = "no", status = 0)
    }
    dispatch_comparisons(comparisons, opts, isar_dir, outdir)
    quit(save = "no", status = 0)
}

features <- features[
    !is.na(features$condition_1) & features$condition_1 == opts$condition_1 &
        !is.na(features$condition_2) & features$condition_2 == opts$condition_2,
    ,
    drop = FALSE
]
if (nrow(features) == 0) {
    stop(
        sprintf("No isoforms found for comparison %s vs %s.", opts$condition_2, opts$condition_1),
        call. = FALSE
    )
}
notes <- c(
    notes,
    sprintf(
        "Comparison: %s (%s vs %s)",
        ifelse(nzchar(opts$comparison_label), opts$comparison_label, sprintf("%s_vs_%s", opts$condition_2, opts$condition_1)),
        opts$condition_2,
        opts$condition_1
    )
)

features$dIF <- safe_num(features$dIF)
features$IF1 <- safe_num(features$IF1)
features$IF2 <- safe_num(features$IF2)
features$isoform_switch_q_value <- safe_num(features$isoform_switch_q_value)
features$gene_switch_q_value <- safe_num(features$gene_switch_q_value)
has_isoform_q_values <- any(!is.na(features$isoform_switch_q_value))
has_gene_q_values <- any(!is.na(features$gene_switch_q_value))
has_statistical_q_values <- has_isoform_q_values || has_gene_q_values
features$neg_log10_q <- safe_neg_log10(features$isoform_switch_q_value)
features$significance <- q_label(features$isoform_switch_q_value, qvalue_cutoff)
features$is_significant_switch <- !is.na(features$isoform_switch_q_value) &
    features$isoform_switch_q_value <= qvalue_cutoff &
    abs(features$dIF) >= dif_cutoff
features$direction <- ifelse(
    is.na(features$dIF) | features$dIF == 0,
    "unchanged",
    ifelse(features$dIF > 0, features$condition_2, features$condition_1)
)
features$isoform_short <- short_isoform_id(features$isoform_id)
features$switch_key <- make_switch_key(features$gene_id, features$condition_1, features$condition_2)
if (!"PTC" %in% colnames(features)) {
    features$PTC <- NA
}

top_switches <- unique(features[, c("gene_id", "gene_name", "condition_1", "condition_2", "gene_switch_q_value")])
if (has_gene_q_values) {
    top_switches <- top_switches[order(top_switches$gene_switch_q_value), , drop = FALSE]
} else {
    gene_rank <- aggregate(list(max_abs_dIF = abs(features$dIF)), by = features[, c("gene_id", "gene_name", "condition_1", "condition_2"), drop = FALSE], FUN = function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
    top_switches <- merge(top_switches, gene_rank, by = c("gene_id", "gene_name", "condition_1", "condition_2"), all.x = TRUE)
    top_switches <- top_switches[order(-top_switches$max_abs_dIF), , drop = FALSE]
}
top_switches$Rank <- seq_len(nrow(top_switches))

if (!"Rank" %in% colnames(top_switches)) {
    top_switches$Rank <- seq_len(nrow(top_switches))
}

top_switches$switch_key <- make_switch_key(top_switches$gene_id, top_switches$condition_1, top_switches$condition_2)
top_switches <- top_switches[!duplicated(top_switches$switch_key), , drop = FALSE]
top_units <- utils::head(top_switches, top_n)
top_keys <- top_units$switch_key
top_genes <- unique(top_units$gene_id)

switch_plot_notes <- character(0)
switch_plots_written <- 0
top_switch_plot_pdf <- file.path(outdir, "top_switch_plots.pdf")
if (!has_statistical_q_values || identical(normalizePath(rds_path, mustWork = FALSE), normalizePath(imported_rds, mustWork = FALSE))) {
    write_text(
        file.path(outdir, "top_switch_plots.txt"),
        "No analyzed switch-test result with statistical q-values was available, so official IsoformSwitchAnalyzeR switchPlot() pages were not written."
    )
} else if (length(top_genes) == 0) {
    write_text(
        file.path(outdir, "top_switch_plots.txt"),
        "No top switch genes were available, so official IsoformSwitchAnalyzeR switchPlot() pages were not written."
    )
} else {
    pdf(top_switch_plot_pdf, width = 12, height = 8)
    cover_labels <- ifelse(!is.na(top_units$gene_name) & nzchar(top_units$gene_name), top_units$gene_name, top_units$gene_id)
    cover_comparisons <- make_comparison_label(top_units$condition_1, top_units$condition_2)
    cover_q <- format_q(top_units$gene_switch_q_value)
    grid::grid.newpage()
    grid::grid.text("Official IsoformSwitchAnalyzeR Top Switch Plots", x = 0.5, y = 0.88, gp = grid::gpar(fontsize = 22, fontface = "bold"))
    grid::grid.text(
        "Selection follows isar_analysis/top_switches.csv from this pipeline run.",
        x = 0.5,
        y = 0.80,
        gp = grid::gpar(fontsize = 13)
    )
    grid::grid.text(
        paste(sprintf("%2d. %s  [%s]  (gene q-value %s)", seq_along(cover_labels), cover_labels, cover_comparisons, cover_q), collapse = "\n"),
        x = 0.08,
        y = 0.66,
        just = c("left", "top"),
        gp = grid::gpar(fontsize = 11, fontfamily = "mono")
    )
    for (i in seq_len(nrow(top_units))) {
        gene <- top_units$gene_id[[i]]
        condition_1 <- top_units$condition_1[[i]]
        condition_2 <- top_units$condition_2[[i]]
        plot_result <- tryCatch({
            IsoformSwitchAnalyzeR::switchPlot(
                switchAnalyzeRlist = switch_list,
                gene = gene,
                condition1 = condition_1,
                condition2 = condition_2
            )
            TRUE
        }, error = function(err) {
            switch_plot_notes <<- c(
                switch_plot_notes,
                sprintf("%s: %s", gene, conditionMessage(err))
            )
            FALSE
        })
        if (isTRUE(plot_result)) {
            switch_plots_written <- switch_plots_written + 1
        }
    }
    dev.off()

    if (switch_plots_written == 0) {
        unlink(top_switch_plot_pdf)
        write_text(
            file.path(outdir, "top_switch_plots.txt"),
            c(
                "IsoformSwitchAnalyzeR switchPlot() could not render any of the selected genes.",
                switch_plot_notes
            )
        )
    } else if (length(switch_plot_notes) > 0) {
        write_text(
            file.path(outdir, "top_switch_plots_warnings.txt"),
            c(
                sprintf("Rendered %d official IsoformSwitchAnalyzeR switchPlot() pages.", switch_plots_written),
                "Some genes could not be rendered:",
                switch_plot_notes
            )
        )
    }
}

top_features <- features[features$switch_key %in% top_keys, , drop = FALSE]
top_features$abs_dIF <- abs(top_features$dIF)
if (has_isoform_q_values) {
    top_features <- top_features[order(top_features$isoform_switch_q_value, -top_features$abs_dIF), , drop = FALSE]
} else {
    top_features <- top_features[order(-top_features$abs_dIF), , drop = FALSE]
}
top_features$dIF_pp <- format_dif_pp(top_features$dIF)
top_features$IF1_percent <- format_if_pct(top_features$IF1)
top_features$IF2_percent <- format_if_pct(top_features$IF2)
top_features$isoform_switch_q_value_formatted <- format_q(top_features$isoform_switch_q_value)
top_features$gene_switch_q_value_formatted <- format_q(top_features$gene_switch_q_value)

candidate_cols <- intersect(
    c("gene_id", "gene_name", "isoform_id", "condition_1", "condition_2", "IF1", "IF2", "dIF", "isoform_switch_q_value", "gene_switch_q_value", "IF1_percent", "IF2_percent", "dIF_pp", "isoform_switch_q_value_formatted", "gene_switch_q_value_formatted", "significance", "direction", "PTC", "iso_biotype"),
    colnames(top_features)
)
write.csv(top_features[, candidate_cols, drop = FALSE], file.path(outdir, "top_isoform_candidates.csv"), row.names = FALSE)

gene_keys <- c("gene_id", "gene_name", "condition_1", "condition_2", "gene_switch_q_value")
gene_isoform_counts <- aggregate(list(isoforms_in_gene = features$isoform_id), by = features[gene_keys], FUN = length)
gene_significant_counts <- aggregate(list(significant_switch_isoforms = features$is_significant_switch), by = features[gene_keys], FUN = function(x) sum(x, na.rm = TRUE))
gene_max_dif <- aggregate(list(max_abs_dIF = abs(features$dIF)), by = features[gene_keys], FUN = function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
gene_summary <- Reduce(function(left, right) merge(left, right, by = gene_keys, all = TRUE), list(gene_isoform_counts, gene_significant_counts, gene_max_dif))
gene_summary$gene_switch_q_value <- safe_num(gene_summary$gene_switch_q_value)
gene_summary$significant_switch_isoforms <- safe_num(gene_summary$significant_switch_isoforms)
gene_summary$max_abs_dIF <- safe_num(gene_summary$max_abs_dIF)
if (has_gene_q_values) {
    gene_summary <- gene_summary[order(gene_summary$gene_switch_q_value, -gene_summary$significant_switch_isoforms, -gene_summary$max_abs_dIF), , drop = FALSE]
} else {
    gene_summary <- gene_summary[order(-gene_summary$significant_switch_isoforms, -gene_summary$max_abs_dIF), , drop = FALSE]
}
gene_summary$gene_switch_q_value_formatted <- format_q(gene_summary$gene_switch_q_value)
gene_summary$max_abs_dIF_pp <- sprintf("%.1f pp", gene_summary$max_abs_dIF * 100)
write.csv(utils::head(gene_summary, top_n), file.path(outdir, "top_gene_summary.csv"), row.names = FALSE)

comparison_pairs <- unique(features[, c("condition_1", "condition_2"), drop = FALSE])
comparison_pairs <- comparison_pairs[stats::complete.cases(comparison_pairs), , drop = FALSE]
if (nrow(comparison_pairs) > 0) {
    volcano_direction_label <- sprintf(
        "Negative dIF: higher isoform usage in %s; positive dIF: higher isoform usage in %s",
        gsub("_", " ", comparison_pairs$condition_1[[1]]),
        gsub("_", " ", comparison_pairs$condition_2[[1]])
    )
} else {
    volcano_direction_label <- "Negative dIF: higher in condition 1; positive dIF: higher in condition 2"
}

top_volcano_labels <- top_features[top_features$is_significant_switch, , drop = FALSE]
top_volcano_labels <- top_volcano_labels[order(top_volcano_labels$isoform_switch_q_value, -abs(top_volcano_labels$dIF)), , drop = FALSE]
top_volcano_labels <- top_volcano_labels[!duplicated(top_volcano_labels$switch_key), , drop = FALSE]
top_volcano_labels <- utils::head(top_volcano_labels, top_n)
top_volcano_labels$label <- ifelse(!is.na(top_volcano_labels$gene_name) & nzchar(top_volcano_labels$gene_name), top_volcano_labels$gene_name, top_volcano_labels$gene_id)

if (has_isoform_q_values) {
    volcano <- ggplot(features, aes(x = dIF, y = neg_log10_q)) +
        geom_point(aes(color = is_significant_switch), alpha = 0.65, size = 1.5) +
        geom_vline(xintercept = c(-dif_cutoff, dif_cutoff), linetype = "dashed", color = "grey55", linewidth = 0.35) +
        geom_hline(yintercept = -log10(qvalue_cutoff), linetype = "dashed", color = "grey55", linewidth = 0.35) +
        geom_text(data = top_volcano_labels, aes(label = label), check_overlap = TRUE, vjust = -0.7, size = 3, color = "grey20") +
        scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey65"), labels = c("TRUE" = cutoff_description, "FALSE" = "other"), name = "Switch candidate") +
        labs(
            title = "Isoform Switch Effect Size vs Statistical Support",
            subtitle = paste(sprintf("Each point is one isoform; dashed lines show %s", cutoff_description), volcano_direction_label, sep = "\n"),
            x = "dIF: change in isoform fraction (condition 2 - condition 1)",
            y = expression(-log[10]("q-value"))
        ) +
        theme_bw(base_size = 12) +
        theme(legend.position = "bottom")
    ggsave(file.path(outdir, "isoform_switch_volcano.png"), volcano, width = 9, height = 6, dpi = 180)
    ggsave(file.path(outdir, "isoform_switch_volcano.pdf"), volcano, width = 9, height = 6)
} else {
    writeLines("No isoform-level q-values were available, so the volcano plot was not written.", file.path(outdir, "isoform_switch_volcano.txt"))
}

top_gene_plot_data <- top_units
top_gene_plot_data$gene_label <- ifelse(!is.na(top_gene_plot_data$gene_name) & nzchar(top_gene_plot_data$gene_name), top_gene_plot_data$gene_name, top_gene_plot_data$gene_id)
top_gene_plot_data$gene_label <- sprintf(
    "%s (%s)",
    top_gene_plot_data$gene_label,
    make_comparison_label(top_gene_plot_data$condition_1, top_gene_plot_data$condition_2)
)
top_gene_plot_data$gene_label <- factor(top_gene_plot_data$gene_label, levels = rev(unique(top_gene_plot_data$gene_label)))
top_gene_plot_data$neg_log10_gene_q <- safe_neg_log10(top_gene_plot_data$gene_switch_q_value)
if (has_gene_q_values) {
    top_gene_plot <- ggplot(top_gene_plot_data, aes(x = neg_log10_gene_q, y = gene_label)) +
        geom_col(fill = "#0072B2", width = 0.75) +
        labs(title = sprintf("Top %d Switching Genes", nrow(top_gene_plot_data)), subtitle = "Ranked by gene-level switch q-value", x = expression(-log[10]("gene switch q-value")), y = NULL) +
        theme_bw(base_size = 12)
    ggsave(file.path(outdir, "top_switching_genes.png"), top_gene_plot, width = 8, height = 6, dpi = 180)
    ggsave(file.path(outdir, "top_switching_genes.pdf"), top_gene_plot, width = 8, height = 6)
} else {
    writeLines("No gene-level q-values were available, so the top switching genes plot was not written.", file.path(outdir, "top_switching_genes.txt"))
}

ptc_data <- features
ptc_data$PTC_status <- ifelse(is.na(ptc_data$PTC), "unknown", ifelse(ptc_data$PTC %in% c(TRUE, "TRUE", "true", "1"), "PTC predicted", "no PTC"))
ptc_data$direction <- ifelse(ptc_data$dIF > 0, paste("higher in", ptc_data$condition_2), ifelse(ptc_data$dIF < 0, paste("higher in", ptc_data$condition_1), "unchanged"))
ptc_data <- ptc_data[ptc_data$is_significant_switch, , drop = FALSE]
if (nrow(ptc_data) > 0) {
    ptc_summary <- as.data.frame(table(ptc_data$PTC_status, ptc_data$direction), stringsAsFactors = FALSE)
    colnames(ptc_summary) <- c("PTC_status", "direction", "n")
    write.csv(ptc_summary, file.path(outdir, "ptc_switch_summary.csv"), row.names = FALSE)
    ptc_plot <- ggplot(ptc_summary, aes(x = direction, y = n, fill = PTC_status)) +
        geom_col(position = "stack") +
        labs(title = "PTC Status Among Significant Isoform Switch Candidates", subtitle = sprintf("Significant means %s", cutoff_description), x = NULL, y = "Number of isoforms", fill = NULL) +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")
    ggsave(file.path(outdir, "ptc_switch_summary.png"), ptc_plot, width = 8, height = 5, dpi = 180)
    ggsave(file.path(outdir, "ptc_switch_summary.pdf"), ptc_plot, width = 8, height = 5)
} else {
    write_empty_csv(file.path(outdir, "ptc_switch_summary.csv"), c("PTC_status", "direction", "n"))
    writeLines(sprintf("No isoforms passed %s, so no PTC summary plot was written.", cutoff_description), file.path(outdir, "ptc_switch_summary.txt"))
}

usage_pdf <- file.path(outdir, "top_gene_isoform_usage.pdf")
usage_dir <- file.path(outdir, "top_gene_isoform_usage")
dir.create(usage_dir, recursive = TRUE, showWarnings = FALSE)
pdf(usage_pdf, width = 10, height = 6)

for (i in seq_len(nrow(top_units))) {
    gene <- top_units$gene_id[[i]]
    switch_key <- top_units$switch_key[[i]]
    gene_df <- features[features$switch_key == switch_key, , drop = FALSE]
    if (nrow(gene_df) == 0) next
    gene_df$max_if <- pmax(gene_df$IF1, gene_df$IF2, na.rm = TRUE)
    if (has_isoform_q_values) {
        gene_df <- gene_df[order(gene_df$isoform_switch_q_value, -abs(gene_df$dIF), -gene_df$max_if), , drop = FALSE]
    } else {
        gene_df <- gene_df[order(-abs(gene_df$dIF), -gene_df$max_if), , drop = FALSE]
    }
    gene_df <- utils::head(gene_df, 10)
    gene_df$isoform_short <- short_isoform_id(gene_df$isoform_id)
    gene_df$annotation_label <- paste0(gene_df$significance, "\nq ", format_q(gene_df$isoform_switch_q_value), "\ndIF ", format_dif_pp(gene_df$dIF))
    gene_label <- gene_df$gene_name[[1]]
    if (is.na(gene_label) || !nzchar(gene_label)) gene_label <- gene
    gene_q_label <- format_q(gene_df$gene_switch_q_value[[1]])
    comparison_label <- make_comparison_label(gene_df$condition_1[[1]], gene_df$condition_2[[1]])
    usage_df <- rbind(
        data.frame(isoform_id = gene_df$isoform_id, isoform_short = gene_df$isoform_short, condition = gene_df$condition_1, IF = gene_df$IF1, stringsAsFactors = FALSE),
        data.frame(isoform_id = gene_df$isoform_id, isoform_short = gene_df$isoform_short, condition = gene_df$condition_2, IF = gene_df$IF2, stringsAsFactors = FALSE)
    )
    usage_df$isoform_short <- factor(usage_df$isoform_short, levels = gene_df$isoform_short[order(gene_df$dIF)])
    label_df <- unique(data.frame(isoform_short = gene_df$isoform_short, annotation_label = gene_df$annotation_label, stringsAsFactors = FALSE))
    max_if <- aggregate(IF ~ isoform_short, usage_df, max, na.rm = TRUE)
    label_df <- merge(label_df, max_if, by = "isoform_short", all.x = TRUE)
    label_df$y <- pmin(label_df$IF + 0.08, 1.05)
    usage_plot <- ggplot(usage_df, aes(x = isoform_short, y = IF, fill = condition)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7) +
        geom_text(data = label_df, aes(x = isoform_short, y = y, label = annotation_label), inherit.aes = FALSE, size = 3.3, lineheight = 0.9) +
        scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.1)) +
        labs(title = sprintf("%s: Isoform Usage", gene_label), subtitle = sprintf("%s; gene switch q-value %s", comparison_label, gene_q_label), x = "Isoform", y = "Isoform fraction", fill = "Condition") +
        theme_bw(base_size = 12) +
        theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
    print(usage_plot)
    usage_prefix <- sprintf("%02d_%s_%s", i, sanitize_filename(gene_label), sanitize_filename(comparison_label))
    ggsave(file.path(usage_dir, sprintf("%s_isoform_usage.png", usage_prefix)), usage_plot, width = 10, height = 6, dpi = 180)
    ggsave(file.path(usage_dir, sprintf("%s_isoform_usage.pdf", usage_prefix)), usage_plot, width = 10, height = 6)
}

dev.off()

writeLines(
    c(
        notes,
        sprintf("RDS file: %s", basename(rds_path)),
        sprintf("Isoforms in feature table: %d", nrow(features)),
        sprintf("Genes in feature table: %d", length(unique(features$gene_id))),
        sprintf("Statistical q-values available: %s", ifelse(has_statistical_q_values, "yes", "no")),
        sprintf("Isoforms passing %s: %d", cutoff_description, sum(features$is_significant_switch, na.rm = TRUE)),
        sprintf("Top gene-comparison units plotted: %d", nrow(top_units)),
        sprintf("Official switchPlot pages written: %d", switch_plots_written)
    ),
    file.path(outdir, "visualization_notes.txt")
)

message("ISAR visualization finished: ", outdir)
