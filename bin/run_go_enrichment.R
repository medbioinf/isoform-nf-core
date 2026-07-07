#!/usr/bin/env Rscript

`%||%` <- function(left, right) {
    if (is.null(left)) right else left
}

parse_args <- function(args) {
    values <- list()
    idx <- 1
    while (idx <= length(args)) {
        key <- args[[idx]]
        if (!startsWith(key, "--")) {
            stop(sprintf("Unexpected positional argument: %s", key), call. = FALSE)
        }
        name <- gsub("-", "_", sub("^--", "", key), fixed = TRUE)
        if (idx == length(args) || startsWith(args[[idx + 1]], "--")) {
            values[[name]] <- TRUE
            idx <- idx + 1
        } else {
            values[[name]] <- args[[idx + 1]]
            idx <- idx + 2
        }
    }
    values
}

required <- function(opts, name) {
    value <- opts[[gsub("-", "_", name, fixed = TRUE)]]
    if (is.null(value) || !nzchar(as.character(value))) {
        stop(sprintf("Missing required argument --%s", name), call. = FALSE)
    }
    as.character(value)
}

clean_ids <- function(ids) {
    ids <- trimws(as.character(ids))
    ids[ids %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
    ids
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

safe_name <- function(x) {
    x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    ifelse(nzchar(x), x, "contrast")
}

choose_column <- function(data, candidates) {
    hits <- candidates[candidates %in% colnames(data)]
    if (length(hits) == 0) NA_character_ else hits[[1]]
}

read_csv_or_empty <- function(path) {
    if (!file.exists(path) || file.info(path)$size == 0) {
        return(data.frame())
    }
    tryCatch(
        read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(err) data.frame()
    )
}

write_empty_enrichment <- function(path) {
    empty <- data.frame(
        contrast = character(),
        geneSet = character(),
        description = character(),
        size = integer(),
        overlap = integer(),
        enrichmentRatio = numeric(),
        pValue = numeric(),
        FDR = numeric(),
        overlapId = character(),
        database = character(),
        stringsAsFactors = FALSE
    )
    write.csv(empty, path, row.names = FALSE)
}

first_non_missing <- function(values) {
    values <- clean_ids(values)
    values <- values[!is.na(values)]
    if (length(values) == 0) NA_character_ else values[[1]]
}

aggregate_gene_scores <- function(rows) {
    if (nrow(rows) == 0) {
        return(data.frame(
            contrast = character(),
            condition_1 = character(),
            condition_2 = character(),
            gene = character(),
            gene_id = character(),
            gene_name = character(),
            min_isoform_switch_q_value = numeric(),
            max_abs_dIF = numeric(),
            significant_isoform_switch = logical(),
            n_isoforms_tested = integer(),
            stringsAsFactors = FALSE
        ))
    }

    keys <- unique(rows[, c("contrast", "condition_1", "condition_2", "gene"), drop = FALSE])
    aggregated <- lapply(seq_len(nrow(keys)), function(idx) {
        key <- keys[idx, , drop = FALSE]
        subset_rows <- rows[
            rows$contrast == key$contrast &
                rows$condition_1 == key$condition_1 &
                rows$condition_2 == key$condition_2 &
                rows$gene == key$gene,
            ,
            drop = FALSE
        ]
        data.frame(
            contrast = key$contrast,
            condition_1 = key$condition_1,
            condition_2 = key$condition_2,
            gene = key$gene,
            gene_id = first_non_missing(subset_rows$gene_id),
            gene_name = first_non_missing(subset_rows$gene_name),
            min_isoform_switch_q_value = min(subset_rows$qvalue, na.rm = TRUE),
            max_abs_dIF = max(abs(subset_rows$dIF), na.rm = TRUE),
            significant_isoform_switch = any(subset_rows$significant, na.rm = TRUE),
            n_isoforms_tested = length(unique(subset_rows$isoform_id)),
            stringsAsFactors = FALSE
        )
    })
    do.call(rbind, aggregated)
}

write_versions <- function(path) {
    r_version <- strsplit(version[["version.string"]], " ")[[1]][3]
    webgestalt_version <- if (requireNamespace("WebGestaltR", quietly = TRUE)) {
        as.character(utils::packageVersion("WebGestaltR"))
    } else {
        "not_available"
    }
    writeLines(
        c(
            '"ISOFORM_GO_ENRICHMENT":',
            sprintf("    r-base: %s", r_version),
            sprintf("    webgestaltr: %s", webgestalt_version)
        ),
        path
    )
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

gene_score_file <- normalizePath(required(opts, "gene-score-file"), mustWork = TRUE)
outdir <- required(opts, "outdir")
organism <- opts[["organism"]] %||% "hsapiens"
gene_id_type <- opts[["gene_id_type"]] %||% "genesymbol"
database <- opts[["database"]] %||% "geneontology_Biological_Process"
qvalue_cutoff <- as.numeric(opts[["qvalue_cutoff"]] %||% 0.05)
dif_cutoff <- as.numeric(opts[["dif_cutoff"]] %||% 0.1)
reference_gene_file <- opts[["reference_gene_file"]]

if (is.na(qvalue_cutoff) || qvalue_cutoff <= 0 || qvalue_cutoff >= 1) {
    stop("--qvalue-cutoff must be between 0 and 1", call. = FALSE)
}
if (is.na(dif_cutoff) || dif_cutoff < 0) {
    stop("--dif-cutoff must be a non-negative number", call. = FALSE)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
go_input_dir <- file.path(outdir, "go_input")
dir.create(go_input_dir, recursive = TRUE, showWarnings = FALSE)
report_root <- file.path(outdir, "webgestalt_report")

notes <- c(
    "GO enrichment summary",
    sprintf("Gene score file: %s", gene_score_file),
    sprintf("Organism: %s", organism),
    sprintf("Gene identifier type: %s", gene_id_type),
    sprintf("Database: %s", database),
    sprintf("q-value cutoff: %.4f", qvalue_cutoff),
    sprintf("dIF cutoff: %.4f", dif_cutoff)
)

scores <- aggregate_gene_scores(data.frame())
combined_enrichment <- data.frame()

scores <- read_csv_or_empty(gene_score_file)
required_cols <- c(
    "contrast",
    "condition_1",
    "condition_2",
    "gene_id",
    "gene_name",
    "min_isoform_switch_q_value",
    "max_abs_dIF",
    "n_isoforms_tested"
)
missing_cols <- setdiff(required_cols, colnames(scores))
if (length(missing_cols) > 0) {
    notes <- c(notes, "Status: skipped", sprintf("Gene score file lacks required columns: %s", paste(missing_cols, collapse = ", ")))
    write.csv(scores, file.path(go_input_dir, "go_gene_scores.csv"), row.names = FALSE)
    write_empty_enrichment(file.path(outdir, "go_enrichment.csv"))
    writeLines(notes, file.path(outdir, "go_summary.txt"))
    write_versions(file.path(outdir, "versions.yml"))
    quit(save = "no", status = 0)
}

scores$gene_id <- clean_ids(scores$gene_id)
scores$gene_name <- clean_ids(scores$gene_name)
use_symbols <- tolower(gene_id_type) %in% c("genesymbol", "gene_symbol", "symbol")
scores$gene <- clean_ids(if (use_symbols) {
    ifelse(!is.na(scores$gene_name), scores$gene_name, scores$gene_id)
} else {
    ifelse(!is.na(scores$gene_id), scores$gene_id, scores$gene_name)
})
scores$min_isoform_switch_q_value <- safe_num(scores$min_isoform_switch_q_value)
scores$max_abs_dIF <- safe_num(scores$max_abs_dIF)
scores$n_isoforms_tested <- as.integer(safe_num(scores$n_isoforms_tested))
scores$significant_isoform_switch <- scores$min_isoform_switch_q_value <= qvalue_cutoff & scores$max_abs_dIF >= dif_cutoff
scores <- scores[!is.na(scores$gene) & !is.na(scores$min_isoform_switch_q_value) & !is.na(scores$max_abs_dIF), , drop = FALSE]
write.csv(scores, file.path(go_input_dir, "go_gene_scores.csv"), row.names = FALSE)

if (nrow(scores) == 0) {
    notes <- c(notes, "Status: skipped", "Gene score file contained no usable gene rows for the configured gene identifier type.")
    write_empty_enrichment(file.path(outdir, "go_enrichment.csv"))
    writeLines(notes, file.path(outdir, "go_summary.txt"))
    write_versions(file.path(outdir, "versions.yml"))
    quit(save = "no", status = 0)
}

if (!is.null(reference_gene_file) && nzchar(as.character(reference_gene_file))) {
    global_reference <- sort(unique(clean_ids(readLines(reference_gene_file, warn = FALSE))))
    global_reference <- global_reference[!is.na(global_reference)]
    notes <- c(notes, sprintf("Using custom GO reference gene file: %s", reference_gene_file))
} else {
    global_reference <- NULL
    notes <- c(notes, "Using all ISAR-tested genes per contrast as the GO reference universe.")
}

run_one_contrast <- function(contrast_name) {
    contrast_scores <- scores[scores$contrast == contrast_name, , drop = FALSE]
    significant_genes <- sort(unique(contrast_scores$gene[contrast_scores$significant_isoform_switch]))
    background_genes <- if (is.null(global_reference)) {
        sort(unique(contrast_scores$gene))
    } else {
        global_reference
    }
    interest_genes <- intersect(significant_genes, background_genes)
    prefix <- safe_name(contrast_name)

    writeLines(significant_genes, file.path(go_input_dir, sprintf("%s_significant_genes.txt", prefix)))
    writeLines(background_genes, file.path(go_input_dir, sprintf("%s_background_genes.txt", prefix)))

    enrichment_path <- file.path(outdir, sprintf("%s_go_enrichment.csv", prefix))
    status <- "completed"
    contrast_notes <- character(0)

    if (length(significant_genes) == 0) {
        status <- "skipped"
        contrast_notes <- c(contrast_notes, "No significant isoform-switching genes passed the configured q-value and dIF cutoffs.")
        write_empty_enrichment(enrichment_path)
    } else if (length(background_genes) < 10) {
        status <- "skipped"
        contrast_notes <- c(contrast_notes, "Skipped WebGestaltR because the GO reference universe has fewer than 10 genes.")
        write_empty_enrichment(enrichment_path)
    } else if (length(interest_genes) == 0) {
        status <- "skipped"
        contrast_notes <- c(contrast_notes, "No significant isoform-switching genes were present in the GO reference universe.")
        write_empty_enrichment(enrichment_path)
    } else {
        result <- tryCatch({
            if (!requireNamespace("WebGestaltR", quietly = TRUE)) {
                stop("The WebGestaltR package is not available in the active runtime")
            }
            report_dir <- file.path(report_root, prefix)
            dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
            WebGestaltR::WebGestaltR(
                enrichMethod = "ORA",
                organism = organism,
                enrichDatabase = database,
                interestGeneFile = file.path(go_input_dir, sprintf("%s_significant_genes.txt", prefix)),
                interestGeneType = gene_id_type,
                referenceGeneFile = file.path(go_input_dir, sprintf("%s_background_genes.txt", prefix)),
                referenceGeneType = gene_id_type,
                minNum = 10,
                maxNum = 500,
                reportNum = 20,
                isOutput = TRUE,
                outputDirectory = report_dir,
                projectName = prefix
            )
        }, error = function(err) {
            status <<- "failed"
            contrast_notes <<- c(contrast_notes, sprintf("WebGestaltR ORA failed: %s", conditionMessage(err)))
            NULL
        })

        if (is.null(result) || nrow(result) == 0) {
            contrast_notes <- c(contrast_notes, "WebGestaltR returned no enriched GO terms.")
            write_empty_enrichment(enrichment_path)
        } else {
            if ("link" %in% colnames(result)) {
                result$link <- NULL
            }
            if (!"database" %in% colnames(result)) {
                result$database <- database
            }
            result$contrast <- contrast_name
            result <- result[, c("contrast", setdiff(colnames(result), "contrast")), drop = FALSE]
            order_columns <- intersect(c("FDR", "pValue"), colnames(result))
            if (length(order_columns) > 0) {
                result <- result[do.call(order, result[order_columns]), , drop = FALSE]
            }
            write.csv(result, enrichment_path, row.names = FALSE)
        }
    }

    enrichment <- read_csv_or_empty(enrichment_path)
    if (nrow(enrichment) > 0 && !"contrast" %in% colnames(enrichment)) {
        enrichment$contrast <- contrast_name
        enrichment <- enrichment[, c("contrast", setdiff(colnames(enrichment), "contrast")), drop = FALSE]
    }

    list(
        contrast = contrast_name,
        status = status,
        significant = length(significant_genes),
        background = length(background_genes),
        interest = length(interest_genes),
        terms = nrow(enrichment),
        notes = contrast_notes,
        enrichment = enrichment
    )
}

contrast_names <- unique(scores$contrast)
results <- lapply(contrast_names, run_one_contrast)
combined_enrichment <- if (length(results) == 0) {
    data.frame()
} else {
    do.call(rbind, lapply(results, function(x) x$enrichment))
}
if (is.null(combined_enrichment) || nrow(combined_enrichment) == 0) {
    write_empty_enrichment(file.path(outdir, "go_enrichment.csv"))
} else {
    write.csv(combined_enrichment, file.path(outdir, "go_enrichment.csv"), row.names = FALSE)
}

summary_lines <- c(
    notes,
    sprintf("Contrasts evaluated: %d", length(results)),
    sprintf("Genes in score table: %d", nrow(scores)),
    sprintf("Enriched terms total: %d", ifelse(is.null(combined_enrichment), 0L, nrow(combined_enrichment)))
)
for (result in results) {
    summary_lines <- c(
        summary_lines,
        sprintf(
            "Contrast %s: status=%s, significant_genes=%d, reference_genes=%d, significant_in_reference=%d, enriched_terms=%d",
            result$contrast,
            result$status,
            result$significant,
            result$background,
            result$interest,
            result$terms
        )
    )
    if (length(result$notes) > 0) {
        summary_lines <- c(summary_lines, paste("- ", unique(result$notes), sep = ""))
    }
}

writeLines(summary_lines, file.path(outdir, "go_summary.txt"))
write_versions(file.path(outdir, "versions.yml"))
