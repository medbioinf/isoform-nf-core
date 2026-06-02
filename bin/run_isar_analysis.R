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
        opts[[sub("^--", "", key)]] <- args[[i + 1]]
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

opts <- parse_args(args)
required <- c("samplesheet", "quant-dir", "gtf", "transcript-fasta", "outdir")
missing <- required[!required %in% names(opts) | !nzchar(unlist(opts[required]))]
if (length(missing) > 0) {
    stop(sprintf("Missing required arguments: %s\n%s", paste(missing, collapse = ", "), usage()), call. = FALSE)
}

samplesheet <- read.csv(opts$samplesheet, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("sample", "condition")
missing_cols <- setdiff(required_cols, colnames(samplesheet))
if (length(missing_cols) > 0) {
    stop(sprintf("Samplesheet is missing required columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

design <- unique(samplesheet[, required_cols, drop = FALSE])
colnames(design) <- c("sampleID", "condition")
if (any(duplicated(design$sampleID))) {
    stop("Each sample must map to exactly one condition for the first ISAR implementation.", call. = FALSE)
}

quant_parent <- normalizePath(opts[["quant-dir"]], mustWork = TRUE)
quant_files <- Sys.glob(file.path(quant_parent, "*", "quant.sf"))
available_dirs <- basename(dirname(quant_files))
design <- design[design$sampleID %in% available_dirs, , drop = FALSE]
if (nrow(design) == 0) {
    stop("No Salmon quantification directories matched samplesheet sample names.", call. = FALSE)
}

dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)
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
comparisons <- NULL
if (nzchar(opts$contrasts)) {
    contrasts <- read.csv(opts$contrasts, stringsAsFactors = FALSE, check.names = FALSE)
    required_contrast_cols <- c("contrast", "case", "control")
    missing_contrast_cols <- setdiff(required_contrast_cols, colnames(contrasts))
    if (length(missing_contrast_cols) > 0) {
        stop(sprintf("Contrast file is missing required columns: %s", paste(missing_contrast_cols, collapse = ", ")), call. = FALSE)
    }
    comparisons <- data.frame(
        condition_1 = contrasts$control,
        condition_2 = contrasts$case,
        stringsAsFactors = FALSE
    )
} else if (length(conditions) == 2) {
    comparisons <- data.frame(condition_1 = conditions[[1]], condition_2 = conditions[[2]], stringsAsFactors = FALSE)
}

message("Creating switchAnalyzeRlist")
switch_list <- importRdata(
    isoformCountMatrix = salmon_quant$counts,
    isoformRepExpression = salmon_quant$abundance,
    designMatrix = design,
    isoformExonAnnoation = filtered_gtf,
    isoformNtFasta = opts[["transcript-fasta"]],
    comparisonsToMake = comparisons,
    ignoreAfterBar = TRUE,
    quiet = TRUE,
    showProgress = FALSE
)

saveRDS(switch_list, file.path(opts$outdir, "switchAnalyzeRlist_imported.rds"))

condition_counts <- table(design$condition)
can_run_test <- length(condition_counts) == 2 && all(condition_counts >= 2) && !is.null(comparisons)
notes <- c(
    sprintf("Samples imported: %d", nrow(design)),
    sprintf("Conditions: %s", paste(names(condition_counts), collapse = ", ")),
    sprintf("Replicates per condition: %s", paste(sprintf("%s=%d", names(condition_counts), as.integer(condition_counts)), collapse = ", "))
)

if (can_run_test) {
    message("Running preFilter() and isoformSwitchTestDEXSeq()")
    switch_list <- preFilter(switch_list)
    result <- tryCatch({
        analyzed <- isoformSwitchTestDEXSeq(
            switchAnalyzeRlist = switch_list,
            alpha = as.numeric(opts$qvalue_cutoff),
            dIFcutoff = as.numeric(opts$dif_cutoff),
            reduceToSwitchingGenes = FALSE,
            quiet = TRUE,
            showProgress = FALSE
        )
        safe_write_csv(extractSwitchSummary(analyzed), file.path(opts$outdir, "switch_summary.csv"))
        safe_write_csv(extractTopSwitches(analyzed, n = as.integer(opts$top_n)), file.path(opts$outdir, "top_switches.csv"))
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
        "Reason: DEXSeq-based testing needs two conditions with at least two samples each."
    )
}

writeLines(notes, file.path(opts$outdir, "analysis_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(opts$outdir, "sessionInfo.txt"))
