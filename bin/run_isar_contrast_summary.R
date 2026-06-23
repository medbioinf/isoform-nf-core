#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
    "Usage:",
    "Rscript run_isar_contrast_summary.R",
    "--isar-dir isar_analysis",
    "--outdir isar_contrast_summary",
    "[--qvalue-cutoff 0.05]",
    "[--dif-cutoff 0.1]",
    "[--top-n 20]",
    sep = "\n"
)

parse_args <- function(args) {
    opts <- list(
        qvalue_cutoff = "0.05",
        dif_cutoff = "0.1",
        top_n = "20"
    )
    i <- 1
    while (i <= length(args)) {
        key <- args[[i]]
        if (!startsWith(key, "--") || i == length(args)) {
            stop(usage, call. = FALSE)
        }
        opts[[gsub("-", "_", sub("^--", "", key), fixed = TRUE)]] <- args[[i + 1]]
        i <- i + 2
    }
    opts
}

write_empty_csv <- function(path, columns) {
    empty <- as.data.frame(setNames(replicate(length(columns), character(0), simplify = FALSE), columns))
    write.csv(empty, path, row.names = FALSE)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

plot_intersections <- function(intersections, matrix_data, contrast_counts, outdir) {
    intersections$intersection_id <- factor(intersections$intersection_id, levels = intersections$intersection_id)
    matrix_data$intersection_id <- factor(matrix_data$intersection_id, levels = levels(intersections$intersection_id))
    contrast_levels <- rev(unique(matrix_data$contrast))
    matrix_data$contrast <- factor(matrix_data$contrast, levels = contrast_levels)
    contrast_counts$contrast <- factor(contrast_counts$contrast, levels = contrast_levels)

    bar_plot <- ggplot(intersections, aes(x = intersection_id, y = n_isoforms)) +
        geom_col(fill = "#595959") +
        geom_text(aes(label = n_isoforms), vjust = -0.35, size = 3.2) +
        labs(title = "Isoform Switch Intersections", y = "Isoform intersections", x = NULL) +
        theme_bw(base_size = 11) +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(face = "bold"),
            plot.margin = margin(5.5, 12, 0, 5.5)
        )

    set_size_plot <- ggplot(contrast_counts, aes(x = n_isoforms, y = contrast)) +
        geom_col(fill = "#595959", width = 0.62) +
        scale_x_reverse(expand = expansion(mult = c(0.05, 0.12))) +
        scale_y_discrete(position = "right", drop = FALSE) +
        labs(x = "Sig. isoforms\nper comparison", y = NULL) +
        theme_bw(base_size = 10) +
        theme(
            panel.grid.major.y = element_blank(),
            axis.text.y = element_text(size = 8),
            axis.ticks.y = element_blank(),
            plot.margin = margin(0, 3, 5.5, 5.5)
        )

    active_data <- matrix_data[matrix_data$active, , drop = FALSE]
    matrix_plot <- ggplot(matrix_data, aes(x = intersection_id, y = contrast)) +
        geom_point(aes(fill = active), shape = 21, size = 3.4, color = "grey75", stroke = 0.4) +
        geom_line(data = active_data, aes(group = intersection_id), linewidth = 0.5, color = "black") +
        scale_fill_manual(values = c("FALSE" = "grey90", "TRUE" = "black"), guide = "none") +
        labs(x = "Intersection rank", y = NULL) +
        theme_bw(base_size = 10) +
        theme(
            panel.grid.major.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            plot.margin = margin(0, 12, 5.5, 5.5)
        )

    pdf(file.path(outdir, "isoform_switch_upset.pdf"), width = 12, height = 8)
    grid::grid.newpage()
    layout <- grid::grid.layout(
        nrow = 2,
        ncol = 2,
        heights = grid::unit(c(0.62, 0.38), "npc"),
        widths = grid::unit(c(0.22, 0.78), "npc")
    )
    grid::pushViewport(grid::viewport(layout = layout))
    print(bar_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(set_size_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(matrix_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
    grid::popViewport()
    dev.off()

    png(file.path(outdir, "isoform_switch_upset.png"), width = 2200, height = 1500, res = 180)
    grid::grid.newpage()
    layout <- grid::grid.layout(
        nrow = 2,
        ncol = 2,
        heights = grid::unit(c(0.62, 0.38), "npc"),
        widths = grid::unit(c(0.22, 0.78), "npc")
    )
    grid::pushViewport(grid::viewport(layout = layout))
    print(bar_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(set_size_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(matrix_plot, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
    grid::popViewport()
    dev.off()
}

opts <- parse_args(args)
required <- c("isar_dir", "outdir")
missing <- required[!required %in% names(opts) | !nzchar(unlist(opts[required]))]
if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s\n%s", paste(missing, collapse = ", "), usage), call. = FALSE)
}

isar_dir <- normalizePath(opts[["isar_dir"]], mustWork = TRUE)
outdir <- opts$outdir
qvalue_cutoff <- as.numeric(opts$qvalue_cutoff)
dif_cutoff <- as.numeric(opts$dif_cutoff)
top_n <- as.integer(opts$top_n)

if (is.na(qvalue_cutoff) || qvalue_cutoff <= 0 || qvalue_cutoff >= 1) {
    stop("--qvalue-cutoff must be between 0 and 1", call. = FALSE)
}
if (is.na(dif_cutoff) || dif_cutoff < 0) {
    stop("--dif-cutoff must be a non-negative number", call. = FALSE)
}
if (is.na(top_n) || top_n < 1) {
    stop("--top-n must be a positive integer", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

analyzed_rds <- file.path(isar_dir, "switchAnalyzeRlist_analyzed.rds")
comparisons_path <- file.path(isar_dir, "comparisons.csv")
notes <- c(
    "ISAR contrast summary",
    sprintf("Input ISAR directory: %s", isar_dir),
    sprintf("q-value cutoff: %.3f", qvalue_cutoff),
    sprintf("dIF cutoff: %.3f", dif_cutoff)
)

if (!file.exists(analyzed_rds)) {
    write_empty_csv(file.path(outdir, "significant_isoform_switches_per_comparison.csv"), c("contrast", "condition_1", "condition_2", "n_isoforms"))
    write_empty_csv(file.path(outdir, "isoform_switch_intersections.csv"), c("intersection_id", "contrasts", "n_isoforms"))
    write_empty_csv(file.path(outdir, "isoform_switch_intersection_members.csv"), c("intersection_id", "contrast", "isoform_id", "gene_id", "gene_name"))
    notes <- c(notes, "No analyzed switchAnalyzeRlist RDS file was found. Contrast summary plots were not generated.")
    writeLines(notes, file.path(outdir, "isar_contrast_summary_notes.txt"))
    quit(save = "no", status = 0)
}

switch_list <- readRDS(analyzed_rds)
features <- switch_list$isoformFeatures
required_cols <- c("isoform_id", "gene_id", "gene_name", "condition_1", "condition_2", "isoform_switch_q_value", "dIF")
missing_cols <- setdiff(required_cols, colnames(features))
if (length(missing_cols) > 0) {
    write_empty_csv(file.path(outdir, "significant_isoform_switches_per_comparison.csv"), c("contrast", "condition_1", "condition_2", "n_isoforms"))
    write_empty_csv(file.path(outdir, "isoform_switch_intersections.csv"), c("intersection_id", "contrasts", "n_isoforms"))
    write_empty_csv(file.path(outdir, "isoform_switch_intersection_members.csv"), c("intersection_id", "contrast", "isoform_id", "gene_id", "gene_name"))
    notes <- c(notes, sprintf("Skipped because isoformFeatures lacks required columns: %s", paste(missing_cols, collapse = ", ")))
    writeLines(notes, file.path(outdir, "isar_contrast_summary_notes.txt"))
    quit(save = "no", status = 0)
}

features$isoform_switch_q_value <- safe_num(features$isoform_switch_q_value)
features$dIF <- safe_num(features$dIF)

contrast_lookup <- unique(features[, c("condition_1", "condition_2"), drop = FALSE])
contrast_lookup$contrast <- sprintf("%s_vs_%s", contrast_lookup$condition_2, contrast_lookup$condition_1)
if (file.exists(comparisons_path)) {
    requested <- read.csv(comparisons_path, stringsAsFactors = FALSE, check.names = FALSE)
    if (all(c("contrast", "condition_1", "condition_2") %in% colnames(requested))) {
        contrast_lookup <- requested[, c("contrast", "condition_1", "condition_2"), drop = FALSE]
    }
}

features <- merge(features, contrast_lookup, by = c("condition_1", "condition_2"), all.x = TRUE)
features$contrast[is.na(features$contrast)] <- sprintf(
    "%s_vs_%s",
    features$condition_2[is.na(features$contrast)],
    features$condition_1[is.na(features$contrast)]
)

significant <- features[
    !is.na(features$isoform_switch_q_value) &
        features$isoform_switch_q_value < qvalue_cutoff &
        !is.na(features$dIF) &
        abs(features$dIF) >= dif_cutoff,
    ,
    drop = FALSE
]
significant <- unique(significant[, c("contrast", "condition_1", "condition_2", "isoform_id", "gene_id", "gene_name"), drop = FALSE])

contrast_order <- contrast_lookup$contrast
contrast_counts <- data.frame(
    contrast = contrast_order,
    stringsAsFactors = FALSE
)
contrast_counts <- merge(
    contrast_counts,
    unique(contrast_lookup[, c("contrast", "condition_1", "condition_2"), drop = FALSE]),
    by = "contrast",
    all.x = TRUE,
    sort = FALSE
)
if (nrow(significant) > 0) {
    counts <- aggregate(isoform_id ~ contrast, significant, function(x) length(unique(x)))
    colnames(counts)[colnames(counts) == "isoform_id"] <- "n_isoforms"
} else {
    counts <- data.frame(contrast = character(), n_isoforms = integer(), stringsAsFactors = FALSE)
}
contrast_counts <- merge(contrast_counts, counts, by = "contrast", all.x = TRUE, sort = FALSE)
contrast_counts$n_isoforms[is.na(contrast_counts$n_isoforms)] <- 0L
contrast_counts <- contrast_counts[match(contrast_order, contrast_counts$contrast), , drop = FALSE]
write.csv(contrast_counts, file.path(outdir, "significant_isoform_switches_per_comparison.csv"), row.names = FALSE)

bar_plot <- ggplot(contrast_counts, aes(x = factor(contrast, levels = contrast), y = n_isoforms)) +
    geom_col(fill = "#4C86B7") +
    geom_text(aes(label = n_isoforms), vjust = -0.3, size = 3.2) +
    labs(
        title = sprintf("Significant isoform switches per comparison (q < %.2g, |dIF| >= %.2g)", qvalue_cutoff, dif_cutoff),
        x = "Comparison",
        y = "Number of significant isoform switches"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(outdir, "significant_isoform_switches_per_comparison.png"), bar_plot, width = 12, height = 6, dpi = 180)
ggsave(file.path(outdir, "significant_isoform_switches_per_comparison.pdf"), bar_plot, width = 12, height = 6)

if (nrow(significant) == 0 || length(unique(significant$contrast)) < 2) {
    write_empty_csv(file.path(outdir, "isoform_switch_intersections.csv"), c("intersection_id", "contrasts", "n_isoforms"))
    write_empty_csv(file.path(outdir, "isoform_switch_intersection_members.csv"), c("intersection_id", "contrast", "isoform_id", "gene_id", "gene_name"))
    notes <- c(
        notes,
        sprintf("Significant isoform rows: %d", nrow(significant)),
        "UpSet plot was not generated because fewer than two comparisons contained significant switches."
    )
    writeLines(notes, file.path(outdir, "isar_contrast_summary_notes.txt"))
    quit(save = "no", status = 0)
}

isoform_meta <- unique(significant[, c("isoform_id", "gene_id", "gene_name"), drop = FALSE])
memberships <- split(significant$contrast, significant$isoform_id)
patterns <- vapply(memberships, function(x) paste(intersect(contrast_order, unique(x)), collapse = "||"), character(1))
intersection_counts <- sort(table(patterns), decreasing = TRUE)
intersection_counts <- utils::head(intersection_counts, top_n)

intersections <- data.frame(
    intersection_id = sprintf("I%02d", seq_along(intersection_counts)),
    contrasts = names(intersection_counts),
    n_isoforms = as.integer(intersection_counts),
    stringsAsFactors = FALSE
)
intersections$contrast_count <- vapply(strsplit(intersections$contrasts, "\\|\\|"), length, integer(1))
intersections <- intersections[order(-intersections$n_isoforms, -intersections$contrast_count, intersections$contrasts), , drop = FALSE]
intersections$intersection_id <- sprintf("I%02d", seq_len(nrow(intersections)))
write.csv(intersections, file.path(outdir, "isoform_switch_intersections.csv"), row.names = FALSE)

member_rows <- do.call(rbind, lapply(seq_len(nrow(intersections)), function(i) {
    pattern <- intersections$contrasts[[i]]
    ids <- names(patterns)[patterns == pattern]
    meta <- isoform_meta[match(ids, isoform_meta$isoform_id), , drop = FALSE]
    data.frame(
        intersection_id = intersections$intersection_id[[i]],
        contrasts = pattern,
        isoform_id = meta$isoform_id,
        gene_id = meta$gene_id,
        gene_name = meta$gene_name,
        stringsAsFactors = FALSE
    )
}))
write.csv(member_rows, file.path(outdir, "isoform_switch_intersection_members.csv"), row.names = FALSE)

matrix_data <- do.call(rbind, lapply(seq_len(nrow(intersections)), function(i) {
    active <- strsplit(intersections$contrasts[[i]], "\\|\\|")[[1]]
    data.frame(
        intersection_id = intersections$intersection_id[[i]],
        contrast = contrast_order,
        active = contrast_order %in% active,
        stringsAsFactors = FALSE
    )
}))
plot_intersections(intersections, matrix_data, contrast_counts, outdir)

notes <- c(
    notes,
    sprintf("Comparisons summarized: %d", nrow(contrast_counts)),
    sprintf("Significant isoform rows: %d", nrow(significant)),
    sprintf("Unique significant isoforms: %d", length(unique(significant$isoform_id))),
    sprintf("Intersections plotted: %d", nrow(intersections))
)
writeLines(notes, file.path(outdir, "isar_contrast_summary_notes.txt"))
