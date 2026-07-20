#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
    paste(
        "Usage: Rscript run_isar_analysis.R",
        "--samplesheet samplesheet.csv",
        "--quant-dir salmon_quant_parent",
        "--gtf annotation.gtf",
        "--transcript-fasta transcripts.fa",
        "--outdir isar_analysis",
        "[--contrasts contrasts.csv]",
        "[--dif-cutoff 0.1]",
        "[--qvalue-cutoff 0.05]",
        "[--top-n 10]"
    )
}

parse_args <- function(args) {
    opts <- list(
        contrasts = "",
        dif_cutoff = "0.1",
        qvalue_cutoff = "0.05",
        top_n = "10"
    )
    i <- 1
    while (i <= length(args)) {
        key <- args[[i]]
        if (!startsWith(key, "--") || i == length(args)) {
            stop(usage(), call. = FALSE)
        }
        option_name <- gsub("-", "_", sub("^--", "", key), fixed = TRUE)
        opts[[option_name]] <- args[[i + 1]]
        i <- i + 2
    }
    opts
}

safe_write_csv <- function(data, path) {
    write.csv(data, path, row.names = FALSE)
}

read_salmon_quantifications <- function(design, quant_parent) {
    count_list <- list()
    abundance_list <- list()

    for (sample_id in design$sampleID) {
        quant_file <- file.path(quant_parent, sample_id, "quant.sf")
        if (!file.exists(quant_file)) {
            stop(sprintf("Missing Salmon quant.sf for sample '%s': %s", sample_id, quant_file), call. = FALSE)
        }
        quant <- read.delim(quant_file, stringsAsFactors = FALSE, check.names = FALSE)
        required_quant_cols <- c("Name", "TPM", "NumReads")
        missing_quant_cols <- setdiff(required_quant_cols, colnames(quant))
        if (length(missing_quant_cols) > 0) {
            stop(
                sprintf("Salmon quant file for '%s' is missing columns: %s", sample_id, paste(missing_quant_cols, collapse = ", ")),
                call. = FALSE
            )
        }
        count_list[[sample_id]] <- stats::setNames(quant$NumReads, quant$Name)
        abundance_list[[sample_id]] <- stats::setNames(quant$TPM, quant$Name)
    }

    isoform_ids <- Reduce(union, lapply(count_list, names))
    counts <- data.frame(isoform_id = isoform_ids, stringsAsFactors = FALSE)
    abundance <- data.frame(isoform_id = isoform_ids, stringsAsFactors = FALSE)
    for (sample_id in design$sampleID) {
        counts[[sample_id]] <- as.numeric(count_list[[sample_id]][isoform_ids])
        abundance[[sample_id]] <- as.numeric(abundance_list[[sample_id]][isoform_ids])
    }
    counts[is.na(counts)] <- 0
    abundance[is.na(abundance)] <- 0

    list(counts = counts, abundance = abundance)
}

extract_gtf_attr <- function(attributes, key) {
    pattern <- sprintf('(?:^|;\\s*)%s "([^"]+)"', key)
    hits <- regmatches(attributes, regexec(pattern, attributes, perl = TRUE))
    vapply(hits, function(match) if (length(match) >= 2) match[[2]] else NA_character_, character(1))
}

normalize_isoform_ids <- function(ids) {
    ids <- as.character(ids)
    ids <- sub("\\s.*$", "", ids)
    ids <- sub("\\|.*$", "", ids)
    ids[!nzchar(ids)] <- NA_character_
    ids
}

write_filtered_gtf <- function(gtf_path, transcript_ids, out_path) {
    transcript_set <- unique(na.omit(normalize_isoform_ids(transcript_ids)))
    if (length(transcript_set) == 0) {
        stop("Cannot filter GTF because no quantified transcript IDs were found.", call. = FALSE)
    }

    gtf <- read.delim(
        gtf_path,
        sep = "\t",
        header = FALSE,
        quote = "",
        comment.char = "#",
        stringsAsFactors = FALSE
    )
    if (ncol(gtf) != 9) {
        stop(sprintf("Expected a 9-column GTF but found %d columns.", ncol(gtf)), call. = FALSE)
    }

    colnames(gtf) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
    gtf$transcript_id <- extract_gtf_attr(gtf$attribute, "transcript_id")
    gtf$gene_id <- extract_gtf_attr(gtf$attribute, "gene_id")

    tx_rows <- !is.na(gtf$transcript_id) & gtf$transcript_id %in% transcript_set
    if (!any(tx_rows)) {
        stop(
            sprintf(
                "No transcript rows in the GTF matched quantified transcript IDs. First IDs: %s",
                paste(utils::head(transcript_set, 5), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    gene_ids <- unique(gtf$gene_id[tx_rows & !is.na(gtf$gene_id)])
    gene_rows <- gtf$feature == "gene" & !is.na(gtf$gene_id) & gtf$gene_id %in% gene_ids
    filtered <- gtf[gene_rows | tx_rows, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")]

    write.table(filtered, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}


safe_num <- function(x) suppressWarnings(as.numeric(x))

clean_ids <- function(ids) {
    ids <- trimws(as.character(ids))
    ids[ids %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
    ids
}

choose_column <- function(data, candidates) {
    hits <- candidates[candidates %in% colnames(data)]
    if (length(hits) == 0) NA_character_ else hits[[1]]
}

first_non_missing <- function(values) {
    values <- clean_ids(values)
    values <- values[!is.na(values)]
    if (length(values) == 0) NA_character_ else values[[1]]
}

normalize_metadata_value <- function(values) {
    if (is.factor(values)) {
        values <- as.character(values)
    }
    if (is.character(values)) {
        values <- trimws(values)
        values[values %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
    }
    values
}

build_design_matrix <- function(samplesheet) {
    sample_ids <- unique(as.character(samplesheet$sample))
    design <- data.frame(sampleID = sample_ids, stringsAsFactors = FALSE)
    design_notes <- character(0)

    collapse_column <- function(column_name) {
        source_values <- normalize_metadata_value(samplesheet[[column_name]])
        collapsed <- lapply(sample_ids, function(sample_id) {
            values <- source_values[as.character(samplesheet$sample) == sample_id]
            values <- unique(values[!is.na(values)])
            if (length(values) > 1) {
                stop(
                    sprintf(
                        "Sample '%s' has conflicting values for design column '%s': %s",
                        sample_id,
                        column_name,
                        paste(values, collapse = ", ")
                    ),
                    call. = FALSE
                )
            }
            if (length(values) == 0) NA else values[[1]]
        })
        if (is.numeric(source_values) || is.integer(source_values)) {
            return(as.numeric(unlist(collapsed)))
        }
        as.character(unlist(collapsed))
    }

    design$condition <- collapse_column("condition")
    if (any(is.na(design$condition))) {
        stop("Every sample must have a non-empty condition.", call. = FALSE)
    }

    technical_columns <- c(
        "sample", "condition", "replicate", "fastq_1", "fastq_2",
        "strandedness", "run_accession", "sra_run", "quant_dir"
    )
    candidate_covariates <- setdiff(colnames(samplesheet), technical_columns)
    if ("batch" %in% candidate_covariates) {
        candidate_covariates <- c("batch", setdiff(candidate_covariates, "batch"))
    }

    included_covariates <- character(0)
    for (column_name in candidate_covariates) {
        values <- collapse_column(column_name)
        present <- !is.na(values)
        if (!any(present)) {
            design_notes <- c(design_notes, sprintf("Design covariate '%s': omitted because it is empty.", column_name))
            next
        }
        if (!all(present)) {
            stop(
                sprintf(
                    "Design covariate '%s' is missing for samples: %s. Covariates must be complete or entirely empty.",
                    column_name,
                    paste(design$sampleID[!present], collapse = ", ")
                ),
                call. = FALSE
            )
        }
        if (length(unique(values)) == 1) {
            design_notes <- c(
                design_notes,
                sprintf("Design covariate '%s': omitted because it is constant across all samples.", column_name)
            )
            next
        }
        if ((is.character(values) || is.factor(values)) && length(unique(values)) == length(values)) {
            stop(
                sprintf(
                    "Design covariate '%s' has a unique categorical value for every sample and would model sample identity rather than a shared effect.",
                    column_name
                ),
                call. = FALSE
            )
        }

        safe_column_name <- make.names(column_name)
        if (safe_column_name %in% colnames(design) || safe_column_name %in% included_covariates) {
            stop(sprintf("Design covariate name '%s' is not unique after R name normalization.", column_name), call. = FALSE)
        }
        design[[safe_column_name]] <- values
        included_covariates <- c(included_covariates, safe_column_name)
        if (!identical(safe_column_name, column_name)) {
            design_notes <- c(
                design_notes,
                sprintf("Design covariate '%s' was normalized to '%s' for R modeling.", column_name, safe_column_name)
            )
        }
    }

    if (length(included_covariates) > 0) {
        model_design <- design
        model_design$condition <- factor(model_design$condition, levels = unique(model_design$condition))
        for (column_name in included_covariates) {
            values <- model_design[[column_name]]
            if (is.numeric(values) || is.integer(values)) {
                if (length(unique(values)) * 2 <= length(values)) {
                    model_design[[column_name]] <- factor(values)
                }
            } else {
                model_design[[column_name]] <- factor(values)
            }
        }
        formula <- stats::reformulate(c("condition", included_covariates), response = NULL, intercept = FALSE)
        model_matrix <- stats::model.matrix(formula, data = model_design)
        if (qr(model_matrix)$rank < ncol(model_matrix)) {
            stop(
                sprintf(
                    paste0(
                        "The design matrix is not full rank after adding covariates: %s. ",
                        "A covariate is likely confounded with condition or another covariate; ",
                        "the effects cannot be estimated independently."
                    ),
                    paste(included_covariates, collapse = ", ")
                ),
                call. = FALSE
            )
        }
        design_notes <- c(
            design_notes,
            sprintf("Design covariates modeled: %s", paste(included_covariates, collapse = ", "))
        )
    } else {
        design_notes <- c(design_notes, "Design covariates modeled: none")
    }

    list(design = design, notes = design_notes)
}

empty_go_gene_scores <- function() {
    data.frame(
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
        qvalue_cutoff = numeric(),
        dif_cutoff = numeric(),
        stringsAsFactors = FALSE
    )
}

write_go_gene_scores <- function(switch_list, contrast_labels, out_path, qvalue_cutoff, dif_cutoff) {
    features <- switch_list$isoformFeatures
    required_cols <- c("condition_1", "condition_2", "isoform_switch_q_value", "dIF")
    if (is.null(features) || nrow(features) == 0 || length(setdiff(required_cols, colnames(features))) > 0) {
        safe_write_csv(empty_go_gene_scores(), out_path)
        return(invisible(FALSE))
    }

    if (!is.null(contrast_labels) && all(c("contrast", "condition_1", "condition_2") %in% colnames(contrast_labels))) {
        contrast_lookup <- unique(contrast_labels[, c("contrast", "condition_1", "condition_2"), drop = FALSE])
    } else {
        contrast_lookup <- unique(features[, c("condition_1", "condition_2"), drop = FALSE])
        contrast_lookup$contrast <- sprintf("%s_vs_%s", contrast_lookup$condition_2, contrast_lookup$condition_1)
        contrast_lookup <- contrast_lookup[, c("contrast", "condition_1", "condition_2"), drop = FALSE]
    }

    features <- merge(features, contrast_lookup, by = c("condition_1", "condition_2"), all.x = TRUE)
    missing_contrast <- is.na(features$contrast)
    features$contrast[missing_contrast] <- sprintf(
        "%s_vs_%s",
        features$condition_2[missing_contrast],
        features$condition_1[missing_contrast]
    )

    gene_id_col <- choose_column(features, c("gene_id", "geneID", "gene"))
    gene_name_col <- choose_column(features, c("gene_name", "geneSymbol", "gene_symbol", "symbol"))
    isoform_id_col <- choose_column(features, c("isoform_id", "transcript_id", "isoform"))

    gene_ids <- if (is.na(gene_id_col)) rep(NA_character_, nrow(features)) else clean_ids(features[[gene_id_col]])
    gene_names <- if (is.na(gene_name_col)) rep(NA_character_, nrow(features)) else clean_ids(features[[gene_name_col]])
    genes <- clean_ids(ifelse(!is.na(gene_names), gene_names, gene_ids))

    rows <- data.frame(
        contrast = features$contrast,
        condition_1 = features$condition_1,
        condition_2 = features$condition_2,
        isoform_id = if (is.na(isoform_id_col)) seq_len(nrow(features)) else features[[isoform_id_col]],
        gene = genes,
        gene_id = gene_ids,
        gene_name = gene_names,
        qvalue = safe_num(features$isoform_switch_q_value),
        dIF = safe_num(features$dIF),
        stringsAsFactors = FALSE
    )
    rows <- rows[!is.na(rows$gene) & !is.na(rows$qvalue) & !is.na(rows$dIF), , drop = FALSE]
    if (nrow(rows) == 0) {
        safe_write_csv(empty_go_gene_scores(), out_path)
        return(invisible(FALSE))
    }
    rows$significant <- rows$qvalue <= qvalue_cutoff & abs(rows$dIF) >= dif_cutoff

    keys <- unique(rows[, c("contrast", "condition_1", "condition_2", "gene"), drop = FALSE])
    scores <- do.call(rbind, lapply(seq_len(nrow(keys)), function(idx) {
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
            qvalue_cutoff = qvalue_cutoff,
            dif_cutoff = dif_cutoff,
            stringsAsFactors = FALSE
        )
    }))
    safe_write_csv(scores, out_path)
    invisible(TRUE)
}

opts <- parse_args(args)
required <- c("samplesheet", "quant_dir", "gtf", "transcript_fasta", "outdir")
missing <- required[!required %in% names(opts) | !nzchar(unlist(opts[required]))]
if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s\n%s", paste(missing, collapse = ", "), usage()), call. = FALSE)
}

qvalue_cutoff <- suppressWarnings(as.numeric(opts$qvalue_cutoff))
dif_cutoff <- suppressWarnings(as.numeric(opts$dif_cutoff))
top_n_numeric <- suppressWarnings(as.numeric(opts$top_n))
if (is.na(qvalue_cutoff) || qvalue_cutoff <= 0 || qvalue_cutoff >= 1) {
    stop("--qvalue-cutoff must be strictly between 0 and 1.", call. = FALSE)
}
if (is.na(dif_cutoff) || dif_cutoff < 0 || dif_cutoff > 1) {
    stop("--dif-cutoff must be between 0 and 1 inclusive.", call. = FALSE)
}
if (is.na(top_n_numeric) || top_n_numeric < 1 || top_n_numeric != floor(top_n_numeric)) {
    stop("--top-n must be a positive integer.", call. = FALSE)
}
top_n <- as.integer(top_n_numeric)

samplesheet <- read.csv(opts$samplesheet, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("sample", "condition")
missing_cols <- setdiff(required_cols, colnames(samplesheet))
if (length(missing_cols) > 0) {
    stop(sprintf("Samplesheet is missing required columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

design_result <- build_design_matrix(samplesheet)
design <- design_result$design
design_notes <- design_result$notes

quant_parent <- normalizePath(opts$quant_dir, mustWork = TRUE)
quant_files <- Sys.glob(file.path(quant_parent, "*", "quant.sf"))
available_dirs <- basename(dirname(quant_files))
design <- design[design$sampleID %in% available_dirs, , drop = FALSE]
if (nrow(design) == 0) {
    stop("No Salmon quantification directories matched samplesheet sample names.", call. = FALSE)
}

dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)
safe_write_csv(empty_go_gene_scores(), file.path(opts$outdir, "go_gene_scores.csv"))
safe_write_csv(design, file.path(opts$outdir, "design_matrix.csv"))
safe_write_csv(
    data.frame(sample = design$sampleID, quant_dir = file.path(quant_parent, design$sampleID), stringsAsFactors = FALSE),
    file.path(opts$outdir, "quant_dirs.csv")
)

message(sprintf("Importing %d Salmon quantifications", nrow(design)))
salmon_quant <- read_salmon_quantifications(design, quant_parent)

quantified_isoforms <- salmon_quant$counts$isoform_id
filtered_gtf <- file.path(opts$outdir, "filtered_annotation.gtf")
write_filtered_gtf(opts$gtf, quantified_isoforms, filtered_gtf)

conditions <- unique(design$condition)
condition_counts <- table(design$condition)
comparisons <- NULL
contrast_labels <- NULL
if (nzchar(opts$contrasts)) {
    contrasts <- read.csv(opts$contrasts, stringsAsFactors = FALSE, check.names = FALSE)
    required_contrast_cols <- c("contrast", "case", "control")
    missing_contrast_cols <- setdiff(required_contrast_cols, colnames(contrasts))
    if (length(missing_contrast_cols) > 0) {
        stop(sprintf("Contrast file is missing required columns: %s", paste(missing_contrast_cols, collapse = ", ")), call. = FALSE)
    }
    missing_conditions <- setdiff(unique(c(contrasts$case, contrasts$control)), conditions)
    if (length(missing_conditions) > 0) {
        stop(sprintf("Contrast file references conditions absent from the samplesheet: %s", paste(missing_conditions, collapse = ", ")), call. = FALSE)
    }
    duplicated_contrasts <- contrasts$contrast[duplicated(contrasts$contrast)]
    if (length(duplicated_contrasts) > 0) {
        stop(sprintf("Contrast names must be unique. Duplicated names: %s", paste(unique(duplicated_contrasts), collapse = ", ")), call. = FALSE)
    }
    comparisons <- data.frame(
        condition_1 = contrasts$control,
        condition_2 = contrasts$case,
        stringsAsFactors = FALSE
    )
    contrast_labels <- data.frame(
        contrast = contrasts$contrast,
        condition_1 = contrasts$control,
        condition_2 = contrasts$case,
        stringsAsFactors = FALSE
    )
} else if (length(conditions) == 2) {
    comparisons <- data.frame(condition_1 = conditions[[1]], condition_2 = conditions[[2]], stringsAsFactors = FALSE)
    contrast_labels <- data.frame(
        contrast = sprintf("%s_vs_%s", comparisons$condition_2, comparisons$condition_1),
        condition_1 = comparisons$condition_1,
        condition_2 = comparisons$condition_2,
        stringsAsFactors = FALSE
    )
}

if (!is.null(contrast_labels)) {
    safe_write_csv(contrast_labels, file.path(opts$outdir, "comparisons.csv"))
}

explicit_covariates <- setdiff(colnames(design), c("sampleID", "condition"))
has_explicit_covariates <- length(explicit_covariates) > 0
import_design <- design[, c("sampleID", "condition"), drop = FALSE]
if (has_explicit_covariates) {
    design_notes <- c(
        design_notes,
        paste0(
            "Explicit covariates are applied in the DEXSeq isoform-usage model. ",
            "IsoformSwitchAnalyzeR automatic surrogate-variable discovery and abundance batch correction were disabled during import."
        )
    )
}

message("Creating switchAnalyzeRlist")
switch_list <- importRdata(
    isoformCountMatrix = salmon_quant$counts,
    isoformRepExpression = salmon_quant$abundance,
    designMatrix = import_design,
    isoformExonAnnoation = filtered_gtf,
    isoformNtFasta = opts$transcript_fasta,
    comparisonsToMake = comparisons,
    detectUnwantedEffects = !has_explicit_covariates,
    ignoreAfterBar = TRUE,
    ignoreAfterPeriod = TRUE,
    quiet = TRUE,
    showProgress = FALSE
)

if (has_explicit_covariates) {
    imported_design <- switch_list$designMatrix
    modeled_design <- design
    modeled_design$condition <- imported_design$condition[
        match(modeled_design$sampleID, imported_design$sampleID)
    ]
    if (any(is.na(modeled_design$condition))) {
        stop("Could not align the validated covariate design with the imported ISAR samples.", call. = FALSE)
    }
    switch_list$designMatrix <- modeled_design
}

saveRDS(switch_list, file.path(opts$outdir, "switchAnalyzeRlist_imported.rds"))

can_run_test <- !is.null(comparisons)
skip_reason <- NULL
if (is.null(comparisons)) {
    skip_reason <- "No explicit contrast file was supplied and the samplesheet does not contain exactly two conditions."
} else {
    contrast_conditions <- unique(c(comparisons$condition_1, comparisons$condition_2))
    low_replicate_conditions <- contrast_conditions[as.integer(condition_counts[contrast_conditions]) < 2]
    if (length(low_replicate_conditions) > 0) {
        can_run_test <- FALSE
        skip_reason <- sprintf(
            "DEXSeq-based testing needs at least two samples per condition used in each contrast. Low-replicate conditions: %s",
            paste(sprintf("%s=%d", low_replicate_conditions, as.integer(condition_counts[low_replicate_conditions])), collapse = ", ")
        )
    }
}
notes <- c(
    sprintf("Samples imported: %d", nrow(design)),
    sprintf("Conditions: %s", paste(names(condition_counts), collapse = ", ")),
    sprintf("Replicates per condition: %s", paste(sprintf("%s=%d", names(condition_counts), as.integer(condition_counts)), collapse = ", ")),
    sprintf("Comparisons requested: %d", ifelse(is.null(comparisons), 0L, nrow(comparisons))),
    design_notes
)

if (can_run_test) {
    message("Running preFilter() and isoformSwitchTestDEXSeq()")
    switch_list <- preFilter(switch_list)
    result <- tryCatch({
        analyzed <- isoformSwitchTestDEXSeq(
            switchAnalyzeRlist = switch_list,
            alpha = qvalue_cutoff,
            dIFcutoff = dif_cutoff,
            reduceToSwitchingGenes = FALSE,
            quiet = TRUE,
            showProgress = FALSE
        )
        safe_write_csv(extractSwitchSummary(analyzed), file.path(opts$outdir, "switch_summary.csv"))
        safe_write_csv(extractTopSwitches(analyzed, n = top_n), file.path(opts$outdir, "top_switches.csv"))
        write_go_gene_scores(
            analyzed,
            contrast_labels,
            file.path(opts$outdir, "go_gene_scores.csv"),
            qvalue_cutoff,
            dif_cutoff
        )
        saveRDS(analyzed, file.path(opts$outdir, "switchAnalyzeRlist_analyzed.rds"))
        list(success = TRUE)
    }, error = function(err) {
        list(success = FALSE, message = conditionMessage(err))
    })

    if (isTRUE(result$success)) {
        notes <- c(notes, "Differential isoform usage test: ran isoformSwitchTestDEXSeq().")
    } else if (grepl("No genes were considered switching", result$message, fixed = TRUE)) {
        notes <- c(notes, "Differential isoform usage test: ran, but no genes passed the current switching cutoffs.")
        writeLines(result$message, file.path(opts$outdir, "no_switches_detected.txt"))
    } else {
        stop(result$message, call. = FALSE)
    }
} else {
    notes <- c(
        notes,
        "Differential isoform usage test: skipped.",
        sprintf("Reason: %s", skip_reason)
    )
}

writeLines(notes, file.path(opts$outdir, "analysis_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(opts$outdir, "sessionInfo.txt"))
