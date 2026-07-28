#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
    "Usage:",
    "Rscript run_annotated_switch_plots.R",
    "<switchAnalyzeRlist_rds_or_dir> <outdir> [n_top] [genes_csv] [condition1] [condition2] [plot_topology] [qvalue_cutoff] [dif_cutoff] [comparisons_csv]",
    "",
    "genes_csv is optional and can be a comma-separated list such as ZNRF3,PBX3.",
    "If genes_csv is '-', the script selects the top n_top genes by gene_switch_q_value.",
    "plot_topology is optional and can be true or false.",
    sep = "\n"
)

if (length(args) < 2 || length(args) > 10) {
    stop(usage, call. = FALSE)
}

first_existing <- function(paths) {
    hits <- paths[file.exists(paths)]
    if (length(hits) == 0) NULL else hits[[1]]
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
rds_path <- if (dir.exists(input_path)) {
    first_existing(file.path(
        input_path,
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
} else {
    input_path
}
if (is.null(rds_path)) {
    stop("No switchAnalyzeRlist RDS file found in: ", input_path, call. = FALSE)
}
outdir <- args[[2]]
n_top <- if (length(args) >= 3) as.integer(args[[3]]) else 10
genes_arg <- if (length(args) >= 4) args[[4]] else "-"
condition1_arg <- if (length(args) >= 5) args[[5]] else NA_character_
condition2_arg <- if (length(args) >= 6) args[[6]] else NA_character_
plot_topology <- if (length(args) >= 7) {
    tolower(args[[7]]) %in% c("true", "t", "1", "yes", "y")
} else {
    TRUE
}
qvalue_cutoff <- if (length(args) >= 8) as.numeric(args[[8]]) else 0.05
dif_cutoff <- if (length(args) >= 9) as.numeric(args[[9]]) else 0.1
comparisons_file <- if (length(args) >= 10 && nzchar(args[[10]])) {
    normalizePath(args[[10]], mustWork = TRUE)
} else {
    ""
}

if (is.na(n_top) || n_top < 1) {
    stop("n_top must be a positive integer.", call. = FALSE)
}
if (is.na(qvalue_cutoff) || qvalue_cutoff <= 0 || qvalue_cutoff >= 1) {
    stop("qvalue_cutoff must be between 0 and 1.", call. = FALSE)
}
if (is.na(dif_cutoff) || dif_cutoff < 0 || dif_cutoff > 1) {
    stop("dif_cutoff must be between 0 and 1 inclusive.", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

safe_filename <- function(x) {
    x <- gsub("[^A-Za-z0-9._-]+", "_", x)
    x[x == ""] <- "unknown"
    x
}

script_path <- function() {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg) != 1) {
        stop("Could not determine the annotated switch plot script path.", call. = FALSE)
    }
    normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}

dispatch_comparisons <- function(features, outdir, comparisons_file) {
    comparisons <- unique(features[, c("condition_1", "condition_2"), drop = FALSE])
    comparisons <- comparisons[
        stats::complete.cases(comparisons) &
            nzchar(comparisons$condition_1) &
            nzchar(comparisons$condition_2),
        ,
        drop = FALSE
    ]
    comparisons$comparison <- sprintf("%s_vs_%s", comparisons$condition_2, comparisons$condition_1)

    if (nzchar(comparisons_file) && file.exists(comparisons_file)) {
        labels <- read.csv(comparisons_file, stringsAsFactors = FALSE, check.names = FALSE)
        required <- c("contrast", "condition_1", "condition_2")
        if (all(required %in% colnames(labels))) {
            for (i in seq_len(nrow(comparisons))) {
                match_row <- labels[
                    labels$condition_1 == comparisons$condition_1[[i]] &
                        labels$condition_2 == comparisons$condition_2[[i]],
                    ,
                    drop = FALSE
                ]
                if (nrow(match_row) > 0 && nzchar(match_row$contrast[[1]])) {
                    comparisons$comparison[[i]] <- match_row$contrast[[1]]
                }
            }
        }
    }
    comparisons$output_directory <- make.unique(safe_filename(comparisons$comparison), sep = "_")

    if (nrow(comparisons) == 0) {
        stop("No tested comparisons were found in the annotated ISAR object.", call. = FALSE)
    }

    statuses <- integer(nrow(comparisons))
    rscript <- file.path(R.home("bin"), "Rscript")
    own_script <- script_path()
    for (i in seq_len(nrow(comparisons))) {
        comparison_outdir <- file.path(outdir, comparisons$output_directory[[i]])
        child_args <- c(
            own_script,
            rds_path,
            comparison_outdir,
            as.character(n_top),
            genes_arg,
            comparisons$condition_1[[i]],
            comparisons$condition_2[[i]],
            ifelse(plot_topology, "true", "false"),
            as.character(qvalue_cutoff),
            as.character(dif_cutoff),
            comparisons_file
        )
        statuses[[i]] <- system2(rscript, args = shQuote(child_args))
    }

    comparisons$status <- ifelse(statuses == 0, "complete", "failed")
    utils::write.csv(
        comparisons,
        file.path(outdir, "annotated_switch_plot_comparisons.csv"),
        row.names = FALSE
    )
    writeLines(
        c(
            "Annotated switch plots by comparison",
            sprintf("Comparisons discovered: %d", nrow(comparisons)),
            sprintf("Comparisons completed: %d", sum(statuses == 0)),
            sprintf("Top genes requested per comparison: %d", n_top),
            if (genes_arg == "-") {
                "Genes were ranked independently within each comparison."
            } else {
                sprintf("Explicit genes requested for every comparison: %s", genes_arg)
            }
        ),
        file.path(outdir, "annotated_switch_plot_notes.txt")
    )

    if (any(statuses != 0)) {
        stop("One or more per-comparison annotated switch plot runs failed. See annotated_switch_plot_comparisons.csv.", call. = FALSE)
    }
}

format_q <- function(q) {
    q <- suppressWarnings(as.numeric(q))
    ifelse(is.na(q), "NA", formatC(q, format = "e", digits = 2))
}

has_nonempty_entry <- function(x, entry_name) {
    entry_name %in% names(x) && !is.null(x[[entry_name]]) && NROW(x[[entry_name]]) > 0
}

has_nonempty_column <- function(x, column_name) {
    column_name %in% colnames(x$isoformFeatures) &&
        any(!is.na(x$isoformFeatures[[column_name]]) & x$isoformFeatures[[column_name]] != "")
}

message("Loading switchAnalyzeRlist: ", rds_path)
switch_list <- readRDS(rds_path)

features <- switch_list$isoformFeatures
required_cols <- c("gene_id", "gene_name", "condition_1", "condition_2", "gene_switch_q_value")
missing_cols <- setdiff(required_cols, colnames(features))
if (length(missing_cols) > 0) {
    stop("isoformFeatures is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

features$gene_switch_q_value <- suppressWarnings(as.numeric(features$gene_switch_q_value))
features$isoform_switch_q_value <- suppressWarnings(as.numeric(features$isoform_switch_q_value))
features$dIF <- suppressWarnings(as.numeric(features$dIF))
features$abs_dIF <- abs(features$dIF)

has_condition1 <- !is.na(condition1_arg) && nzchar(condition1_arg)
has_condition2 <- !is.na(condition2_arg) && nzchar(condition2_arg)
if (xor(has_condition1, has_condition2)) {
    stop("condition1 and condition2 must be supplied together.", call. = FALSE)
}
if (!has_condition1) {
    dispatch_comparisons(features, outdir, comparisons_file)
    quit(save = "no", status = 0)
}

condition1 <- if (!is.na(condition1_arg) && nzchar(condition1_arg)) condition1_arg else features$condition_1[[1]]
condition2 <- if (!is.na(condition2_arg) && nzchar(condition2_arg)) condition2_arg else features$condition_2[[1]]

matches_comparison <-
    !is.na(features$condition_1) & features$condition_1 == condition1 &
    !is.na(features$condition_2) & features$condition_2 == condition2
comparison_features <- features[matches_comparison, , drop = FALSE]
if (nrow(comparison_features) == 0) {
    stop(sprintf("No isoforms found for comparison %s vs %s.", condition1, condition2), call. = FALSE)
}

if (!is.na(genes_arg) && nzchar(genes_arg) && genes_arg != "-") {
    genes_to_plot <- trimws(strsplit(genes_arg, ",", fixed = TRUE)[[1]])
    genes_to_plot <- genes_to_plot[nzchar(genes_to_plot)]
} else {
    significant_features <- comparison_features[
        !is.na(comparison_features$isoform_switch_q_value) &
            comparison_features$isoform_switch_q_value <= qvalue_cutoff &
            !is.na(comparison_features$abs_dIF) &
            comparison_features$abs_dIF >= dif_cutoff,
        ,
        drop = FALSE
    ]
    gene_summary <- unique(significant_features[, c("gene_id", "gene_name", "gene_switch_q_value"), drop = FALSE])
    has_gene_identifier <-
        (!is.na(gene_summary$gene_id) & nzchar(gene_summary$gene_id)) |
        (!is.na(gene_summary$gene_name) & nzchar(gene_summary$gene_name))
    gene_summary <- gene_summary[has_gene_identifier & !is.na(gene_summary$gene_switch_q_value), , drop = FALSE]
    gene_summary <- gene_summary[order(gene_summary$gene_switch_q_value), , drop = FALSE]
    gene_summary$plot_identifier <- ifelse(
        !is.na(gene_summary$gene_id) & nzchar(gene_summary$gene_id),
        gene_summary$gene_id,
        gene_summary$gene_name
    )
    gene_summary <- gene_summary[!duplicated(gene_summary$plot_identifier), , drop = FALSE]
    genes_to_plot <- head(gene_summary$plot_identifier, n_top)
}

gene_summary <- unique(comparison_features[, c("gene_id", "gene_name", "gene_switch_q_value"), drop = FALSE])
gene_summary <- gene_summary[gene_summary$gene_id %in% genes_to_plot | gene_summary$gene_name %in% genes_to_plot, , drop = FALSE]
gene_summary <- gene_summary[order(gene_summary$gene_switch_q_value), , drop = FALSE]
gene_summary$gene_switch_q_value_formatted <- format_q(gene_summary$gene_switch_q_value)
utils::write.csv(gene_summary, file.path(outdir, "annotated_switch_plot_genes.csv"), row.names = FALSE)

annotation_status <- data.frame(
    annotation = c(
        "ORF/PTC",
        "Pfam protein domains",
        "SignalP signal peptide",
        "IUPred2A/NetSurfP IDR",
        "DeepLoc2 subcellular location",
        "DeepTMHMM topology"
    ),
    status = c(
        ifelse(has_nonempty_entry(switch_list, "orfAnalysis"), "available", "missing"),
        ifelse(has_nonempty_entry(switch_list, "domainAnalysis"), "available", "missing"),
        ifelse(has_nonempty_entry(switch_list, "signalPeptideAnalysis"), "available", "missing"),
        ifelse(has_nonempty_entry(switch_list, "idrAnalysis"), "available", "missing"),
        ifelse(has_nonempty_column(switch_list, "sub_cell_location"), "available", "missing"),
        ifelse(has_nonempty_entry(switch_list, "topologyAnalysis"), "available", "missing")
    ),
    stringsAsFactors = FALSE
)
utils::write.csv(annotation_status, file.path(outdir, "annotation_status.csv"), row.names = FALSE)

plot_records <- list()
for (i in seq_along(genes_to_plot)) {
    gene <- genes_to_plot[[i]]
    matches_gene <-
        (!is.na(comparison_features$gene_id) & comparison_features$gene_id == gene) |
        (!is.na(comparison_features$gene_name) & comparison_features$gene_name == gene)
    gene_rows <- comparison_features[matches_gene, , drop = FALSE]
    if (nrow(gene_rows) == 0) {
        warning("Skipping gene not found in comparison: ", gene)
        next
    }

    gene_id_candidates <- gene_rows$gene_id[!is.na(gene_rows$gene_id) & nzchar(gene_rows$gene_id)]
    gene_name_candidates <- gene_rows$gene_name[!is.na(gene_rows$gene_name) & nzchar(gene_rows$gene_name)]
    gene_id <- if (length(gene_id_candidates) > 0) gene_id_candidates[[1]] else gene
    gene_name <- if (length(gene_name_candidates) > 0) gene_name_candidates[[1]] else NA_character_
    label <- if (!is.na(gene_name) && nzchar(gene_name)) gene_name else gene_id
    file_prefix <- sprintf("%02d_%s_annotated_switch", i, safe_filename(label))
    pdf_path <- file.path(outdir, paste0(file_prefix, ".pdf"))
    png_path <- file.path(outdir, paste0(file_prefix, ".png"))

    draw_plot <- function(use_topology) {
        switchPlot(
            switchAnalyzeRlist = switch_list,
            gene = gene_id,
            condition1 = condition1,
            condition2 = condition2,
            IFcutoff = 0.05,
            dIFcutoff = dif_cutoff,
            alphas = c(qvalue_cutoff, min(qvalue_cutoff, 0.001)),
            plotTopology = use_topology,
            localTheme = ggplot2::theme_bw(base_size = 12)
        )
    }

    render_plot <- function(device_label) {
        tryCatch({
            draw_plot(plot_topology)
        }, error = function(err) {
            if (plot_topology) {
                warning(
                    "Could not plot ", device_label, " for ", label,
                    " with topology enabled. Retrying without topology: ",
                    conditionMessage(err)
                )
                tryCatch({
                    draw_plot(FALSE)
                }, error = function(fallback_err) {
                    plot.new()
                    title(main = paste("Could not plot", label))
                    text(0.5, 0.5, conditionMessage(fallback_err))
                    warning("Could not plot ", device_label, " for ", label, ": ", conditionMessage(fallback_err))
                })
            } else {
                plot.new()
                title(main = paste("Could not plot", label))
                text(0.5, 0.5, conditionMessage(err))
                warning("Could not plot ", device_label, " for ", label, ": ", conditionMessage(err))
            }
        })
    }

    message("Plotting ", label, " -> ", pdf_path)
    pdf(file = pdf_path, onefile = FALSE, width = 12, height = 7)
    render_plot("PDF")
    dev.off()

    png(filename = png_path, width = 2400, height = 1400, res = 200)
    render_plot("PNG")
    dev.off()

    best_isoform <- gene_rows[order(gene_rows$isoform_switch_q_value, -gene_rows$abs_dIF), , drop = FALSE][1, ]
    plot_records[[length(plot_records) + 1]] <- data.frame(
        gene_id = gene_id,
        gene_name = label,
        condition_1 = condition1,
        condition_2 = condition2,
        gene_switch_q_value = gene_rows$gene_switch_q_value[[1]],
        gene_switch_q_value_formatted = format_q(gene_rows$gene_switch_q_value[[1]]),
        strongest_isoform_id = best_isoform$isoform_id,
        strongest_isoform_q_value = best_isoform$isoform_switch_q_value,
        strongest_isoform_q_value_formatted = format_q(best_isoform$isoform_switch_q_value),
        strongest_isoform_dIF = best_isoform$dIF,
        pdf_file = basename(pdf_path),
        png_file = basename(png_path),
        stringsAsFactors = FALSE
    )
}

plot_table <- if (length(plot_records) > 0) do.call(rbind, plot_records) else NULL
if (!is.null(plot_table) && nrow(plot_table) > 0) {
    utils::write.csv(plot_table, file.path(outdir, "annotated_switch_plot_summary.csv"), row.names = FALSE)
} else {
    utils::write.csv(
        data.frame(
            gene_id = character(),
            gene_name = character(),
            condition_1 = character(),
            condition_2 = character(),
            gene_switch_q_value = numeric(),
            pdf_file = character(),
            png_file = character()
        ),
        file.path(outdir, "annotated_switch_plot_summary.csv"),
        row.names = FALSE
    )
}

notes <- c(
    "Annotated ISAR switch plot summary",
    sprintf("Input RDS: %s", rds_path),
    sprintf("Output directory: %s", normalizePath(outdir, mustWork = FALSE)),
    sprintf("Comparison: %s vs %s", condition1, condition2),
    sprintf("Switch q-value cutoff: %s", format(qvalue_cutoff)),
    sprintf("Switch dIF cutoff: %s", format(dif_cutoff)),
    sprintf("Topology plotting requested: %s", plot_topology),
    sprintf("Genes requested/plotted: %s", paste(genes_to_plot, collapse = ", ")),
    "",
    "Annotation layers:",
    sprintf("- %s: %s", annotation_status$annotation, annotation_status$status),
    "",
    "Interpretation:",
    "These plots are generated with IsoformSwitchAnalyzeR::switchPlot().",
    "Protein domains are shown when analyzePFAM() results are present.",
    "Topology requires analyzeDeepTMHMM() results. Subcellular location requires analyzeDeepLoc2() results.",
    "IDR tracks require analyzeIUPred2A() or analyzeNetSurfP2() results. Signal peptides require analyzeSignalP() results."
)
writeLines(notes, file.path(outdir, "annotated_switch_plot_notes.txt"))

message("Wrote annotated switch plots to: ", outdir)
