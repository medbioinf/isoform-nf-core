#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
    "Usage:",
    "Rscript run_pfam_visualization.R",
    "--pfam-dir pfam_import",
    "--outdir pfam_visualization",
    "[--top-n 12]",
    sep = "\n"
)

parse_args <- function(args) {
    opts <- list(top_n = "12")
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

write_empty_csv <- function(path, columns) {
    empty <- as.data.frame(setNames(replicate(length(columns), character(0), simplify = FALSE), columns))
    write.csv(empty, path, row.names = FALSE)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

format_q <- function(q) {
    q <- safe_num(q)
    ifelse(is.na(q), "NA", formatC(q, format = "e", digits = 2))
}

format_condition <- function(condition) {
    gsub("_", " ", as.character(condition))
}

first_existing_column <- function(data, candidates) {
    for (candidate in candidates) {
        if (candidate %in% colnames(data)) {
            return(data[[candidate]])
        }
    }
    rep(NA_character_, nrow(data))
}

domain_label <- function(values) {
    values <- unique(values[!is.na(values) & values != ""])
    if (length(values) == 0) {
        return("none")
    }
    paste(values, collapse = "; ")
}

opts <- parse_args(args)
required <- c("pfam-dir", "outdir")
missing <- required[!required %in% names(opts) | !nzchar(unlist(opts[required]))]
if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s\n%s", paste(missing, collapse = ", "), usage()), call. = FALSE)
}

pfam_dir <- normalizePath(opts[["pfam-dir"]], mustWork = TRUE)
outdir <- opts$outdir
top_n <- as.integer(opts$top_n)
if (is.na(top_n) || top_n < 1) {
    stop("--top-n must be a positive integer", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
domain_plot_dir <- file.path(outdir, "domain_architecture")
dir.create(domain_plot_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- file.path(
    pfam_dir,
    c("pfam_domain_analysis.csv", "pfam_switch_consequences.csv", "significant_switches_with_pfam.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
    write_empty_csv(file.path(outdir, "pfam_top_domain_change_candidates.csv"), c("gene_id", "gene_name", "switchConsequence"))
    write_empty_csv(file.path(outdir, "pfam_switch_consequence_summary.csv"), c("switch_consequence", "usage_shift", "n_switches"))
    writeLines(
        c(
            "Pfam visualization summary",
            sprintf("Input Pfam import directory: %s", pfam_dir),
            sprintf("Skipped because required files are missing: %s", paste(basename(missing_files), collapse = ", "))
        ),
        file.path(outdir, "pfam_visualization_notes.txt")
    )
    quit(save = "no", status = 0)
}

domains <- read.csv(file.path(pfam_dir, "pfam_domain_analysis.csv"), stringsAsFactors = FALSE)
consequences <- read.csv(file.path(pfam_dir, "pfam_switch_consequences.csv"), stringsAsFactors = FALSE)
switches <- read.csv(file.path(pfam_dir, "significant_switches_with_pfam.csv"), stringsAsFactors = FALSE)

domains$orf_aa_start <- safe_num(domains$orf_aa_start)
domains$orf_aa_end <- safe_num(domains$orf_aa_end)
domains$significant <- safe_num(domains$significant)
domains <- domains[is.na(domains$significant) | domains$significant == 1, , drop = FALSE]

switches$isoform_switch_q_value <- safe_num(switches$isoform_switch_q_value)
switches$dIF <- safe_num(switches$dIF)
switches$abs_dIF <- abs(switches$dIF)

domains_for_isoform <- function(isoform_id) {
    rows <- domains[domains$isoform_id == isoform_id, , drop = FALSE]
    rows <- rows[order(rows$orf_aa_start, rows$orf_aa_end, rows$hmm_name), , drop = FALSE]
    rows
}

domain_names_for_isoform <- function(isoform_id) {
    domain_label(domains_for_isoform(isoform_id)$hmm_name)
}

consequences$switchConsequence <- trimws(as.character(consequences$switchConsequence))
consequences$switchConsequence[consequences$switchConsequence %in% c("", "NA", "NaN")] <- NA

domain_changes <- consequences[!is.na(consequences$switchConsequence), , drop = FALSE]
domain_changes <- domain_changes[order(domain_changes$switchConsequence, domain_changes$gene_name), , drop = FALSE]

if (nrow(domain_changes) == 0) {
    write_empty_csv(file.path(outdir, "pfam_top_domain_change_candidates.csv"), c("gene_id", "gene_name", "switchConsequence"))
    write_empty_csv(file.path(outdir, "pfam_switch_consequence_summary.csv"), c("switch_consequence", "usage_shift", "n_switches"))
    writeLines(
        c(
            "Pfam visualization summary",
            sprintf("Input Pfam import directory: %s", pfam_dir),
            "No Pfam domain gain/loss/switch consequences were found."
        ),
        file.path(outdir, "pfam_visualization_notes.txt")
    )
    quit(save = "no", status = 0)
}

best_switch_by_gene <- switches[order(switches$isoform_switch_q_value, -switches$abs_dIF), , drop = FALSE]
best_switch_by_gene <- best_switch_by_gene[!duplicated(best_switch_by_gene$gene_id), , drop = FALSE]
best_switch_cols <- intersect(
    c("gene_id", "condition_1", "condition_2", "isoform_switch_q_value", "gene_switch_q_value", "dIF", "dIF_pp"),
    colnames(best_switch_by_gene)
)
best_switch_by_gene <- best_switch_by_gene[, best_switch_cols, drop = FALSE]

candidate_table <- merge(domain_changes, best_switch_by_gene, by = "gene_id", all.x = TRUE)
candidate_table$q_value_formatted <- format_q(candidate_table$isoform_switch_q_value)
candidate_table$condition_1_effect <- first_existing_column(candidate_table, c("condition_1.x", "condition_1.y", "condition_1"))
candidate_table$condition_2_effect <- first_existing_column(candidate_table, c("condition_2.x", "condition_2.y", "condition_2"))
candidate_table$upregulated_condition <- ifelse(
    is.na(candidate_table$dIF),
    NA_character_,
    ifelse(candidate_table$dIF >= 0, candidate_table$condition_2_effect, candidate_table$condition_1_effect)
)
candidate_table$downregulated_condition <- ifelse(
    is.na(candidate_table$dIF),
    NA_character_,
    ifelse(candidate_table$dIF >= 0, candidate_table$condition_1_effect, candidate_table$condition_2_effect)
)
candidate_table$usage_shift <- ifelse(
    is.na(candidate_table$upregulated_condition) | is.na(candidate_table$downregulated_condition),
    "unknown direction",
    paste(format_condition(candidate_table$downregulated_condition), "to", format_condition(candidate_table$upregulated_condition))
)
candidate_table$upregulated_domains <- vapply(candidate_table$isoformUpregulated, domain_names_for_isoform, character(1))
candidate_table$downregulated_domains <- vapply(candidate_table$isoformDownregulated, domain_names_for_isoform, character(1))

candidate_table <- candidate_table[order(
    candidate_table$isoform_switch_q_value,
    -abs(safe_num(candidate_table$dIF)),
    candidate_table$gene_name
), , drop = FALSE]

utils::write.csv(
    candidate_table,
    file.path(outdir, "pfam_top_domain_change_candidates.csv"),
    row.names = FALSE
)

summary_table <- aggregate(
    list(n_switches = candidate_table$switchConsequence),
    by = candidate_table[, c("switchConsequence", "usage_shift"), drop = FALSE],
    FUN = length
)
colnames(summary_table)[colnames(summary_table) == "switchConsequence"] <- "switch_consequence"
summary_table <- summary_table[order(-summary_table$n_switches, summary_table$switch_consequence, summary_table$usage_shift), , drop = FALSE]
utils::write.csv(summary_table, file.path(outdir, "pfam_switch_consequence_summary.csv"), row.names = FALSE)

summary_plot <- ggplot(summary_table, aes(x = reorder(switch_consequence, n_switches, FUN = sum), y = n_switches, fill = usage_shift)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = n_switches), position = position_stack(vjust = 0.5), size = 3.6, color = "white") +
    coord_flip(clip = "off") +
    labs(
        title = "Pfam Domain Consequences",
        subtitle = "Significant isoform switches where protein-domain content changes; fill shows the isoform usage shift",
        x = NULL,
        y = "Number of switch events",
        fill = "Usage shift"
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title.position = "plot",
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(8, 24, 8, 8)
    )

ggsave(file.path(outdir, "pfam_switch_consequence_summary.png"), summary_plot, width = 8, height = 4.5, dpi = 180)
ggsave(file.path(outdir, "pfam_switch_consequence_summary.pdf"), summary_plot, width = 8, height = 4.5)

plot_domain_architecture <- function(row) {
    up_id <- row[["isoformUpregulated"]]
    down_id <- row[["isoformDownregulated"]]
    gene_name <- row[["gene_name"]]
    consequence <- row[["switchConsequence"]]
    q_label <- row[["q_value_formatted"]]
    dif_label <- if ("dIF_pp" %in% names(row) && !is.na(row[["dIF_pp"]])) row[["dIF_pp"]] else ""
    up_condition <- if ("upregulated_condition" %in% names(row)) row[["upregulated_condition"]] else NA_character_
    down_condition <- if ("downregulated_condition" %in% names(row)) row[["downregulated_condition"]] else NA_character_
    usage_shift <- if ("usage_shift" %in% names(row)) row[["usage_shift"]] else "unknown direction"

    roles <- c(
        paste0("Downregulated isoform", ifelse(is.na(down_condition), "", paste0(" (higher in ", format_condition(down_condition), ")"))),
        paste0("Upregulated isoform", ifelse(is.na(up_condition), "", paste0(" (higher in ", format_condition(up_condition), ")")))
    )

    up_domains <- domains_for_isoform(up_id)
    down_domains <- domains_for_isoform(down_id)
    up_domains$isoform_role <- character(nrow(up_domains))
    down_domains$isoform_role <- character(nrow(down_domains))
    if (nrow(up_domains) > 0) {
        up_domains$isoform_role <- roles[[2]]
    }
    if (nrow(down_domains) > 0) {
        down_domains$isoform_role <- roles[[1]]
    }
    plot_domains <- rbind(up_domains, down_domains)
    max_end <- suppressWarnings(max(plot_domains$orf_aa_end, na.rm = TRUE))
    if (!is.finite(max_end)) {
        max_end <- 100
    }

    base <- data.frame(
        isoform_role = roles,
        y = seq_along(roles),
        x_start = 1,
        x_end = max_end,
        isoform_id = c(down_id, up_id),
        stringsAsFactors = FALSE
    )

    plot_domains$y <- match(plot_domains$isoform_role, roles)
    plot_domains$label <- ifelse(
        is.na(plot_domains$hmm_name) | plot_domains$hmm_name == "",
        plot_domains$hmm_acc,
        plot_domains$hmm_name
    )

    p <- ggplot() +
        geom_segment(
            data = base,
            aes(x = x_start, xend = x_end, y = y, yend = y),
            linewidth = 3,
            color = "#d7d1c8",
            lineend = "round"
        ) +
        geom_rect(
            data = plot_domains,
            aes(
                xmin = pmax(orf_aa_start, 1),
                xmax = pmax(orf_aa_end, orf_aa_start + 1),
                ymin = y - 0.18,
                ymax = y + 0.18,
                fill = label
            ),
            color = "white",
            linewidth = 0.3
        ) +
        geom_text(
            data = plot_domains,
            aes(
                x = (orf_aa_start + orf_aa_end) / 2,
                y = y,
                label = label
            ),
            size = 2.8,
            color = "white",
            check_overlap = TRUE
        ) +
        geom_text(
            data = base,
            aes(x = 1, y = y, label = isoform_id),
            hjust = 0,
            nudge_y = 0.32,
            size = 3.1,
            color = "#4b4b4b"
        ) +
        scale_y_continuous(breaks = seq_along(roles), labels = roles, limits = c(0.55, length(roles) + 0.55)) +
        labs(
            title = paste0(gene_name, ": ", consequence),
            subtitle = paste0("Usage shift: ", usage_shift, "; q-value ", q_label, ifelse(dif_label == "", "", paste0(", dIF ", dif_label))),
            x = "Amino-acid position in predicted ORF",
            y = NULL,
            fill = "Pfam domain"
        ) +
        theme_bw(base_size = 11) +
        theme(
            plot.title.position = "plot",
            legend.position = "bottom",
            panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank()
        )

    if (nrow(plot_domains) == 0) {
        p <- p + annotate(
            "text",
            x = max_end / 2,
            y = 1.5,
            label = "No significant Pfam domains for these isoforms",
            color = "#555555"
        )
    }

    p
}

top_candidates <- head(candidate_table, top_n)
safe_gene <- gsub("[^A-Za-z0-9._-]+", "_", top_candidates$gene_name)
safe_gene[safe_gene == ""] <- paste0("candidate_", seq_len(sum(safe_gene == "")))
top_candidates$plot_file_prefix <- sprintf("%02d_%s", seq_len(nrow(top_candidates)), safe_gene)

pdf(file.path(outdir, "pfam_top_domain_architectures.pdf"), width = 10, height = 5.6)
for (i in seq_len(nrow(top_candidates))) {
    p <- plot_domain_architecture(top_candidates[i, , drop = FALSE])
    print(p)
    ggsave(
        file.path(domain_plot_dir, paste0(top_candidates$plot_file_prefix[[i]], "_domain_architecture.png")),
        p,
        width = 10,
        height = 5.6,
        dpi = 180
    )
    ggsave(
        file.path(domain_plot_dir, paste0(top_candidates$plot_file_prefix[[i]], "_domain_architecture.pdf")),
        p,
        width = 10,
        height = 5.6
    )
}
dev.off()

notes <- c(
    "Pfam visualization summary",
    sprintf("Input Pfam import directory: %s", pfam_dir),
    sprintf("Top candidates plotted: %d", nrow(top_candidates)),
    sprintf("Domain consequence events: %d", nrow(domain_changes)),
    sprintf(
        "Consequence counts: %s",
        paste(sprintf("%s (%s)=%d", summary_table$switch_consequence, summary_table$usage_shift, summary_table$n_switches), collapse = ", ")
    ),
    "",
    "Interpretation:",
    "Domain gain/loss/switch means that the isoforms involved in a significant switch differ in predicted Pfam protein-domain content.",
    "The architecture plots compare the upregulated and downregulated isoform for a gene. Colored blocks are Pfam domains along the predicted protein sequence.",
    "The usage-shift labels describe which condition loses usage of the downregulated isoform and which condition gains usage of the upregulated isoform."
)
writeLines(notes, file.path(outdir, "pfam_visualization_notes.txt"))

message("Wrote Pfam visualization outputs to: ", outdir)
