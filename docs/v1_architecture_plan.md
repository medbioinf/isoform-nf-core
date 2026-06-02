# V1 Architecture Status

Last updated: 2026-06-02

This document summarizes the current V1 architecture after the reusable workflow path was implemented.

## Current Workflow

The implemented workflow is:

`FASTQ or SRA manifest -> FastQC -> cat/fastq -> fastp -> Salmon index / quant -> IsoformSwitchAnalyzeR -> MultiQC`

The top-level input branch supports either a local FASTQ samplesheet or an SRA manifest. After SRA runs are downloaded and converted to FASTQ, both input modes continue through the same downstream path.

## Reused nf-core Modules

The pipeline currently reuses these nf-core modules:

- `fastqc`
- `cat/fastq`
- `fastp`
- `salmon/index`
- `salmon/quant`
- `multiqc`

This keeps the standard RNA-seq processing steps close to nf-core conventions instead of reimplementing them locally.

## Local Modules

### `FETCH_SRA_FASTQ`

Purpose:

- download SRA/ENA/DRA runs with `prefetch`
- convert archives with `fasterq-dump`
- compress FASTQs with `pigz`
- normalize output names for downstream processing

Reason for local implementation:

- the project needs a small manifest-driven SRA mode that feeds directly into the same internal FASTQ path
- the first real-data tests showed that direct `fasterq-dump SRR...` was less robust than `prefetch` followed by local conversion

### `ISAR_ANALYSIS`

Purpose:

- collect Salmon quantification directories
- build count and TPM matrices
- filter the GTF to quantified transcripts
- import into IsoformSwitchAnalyzeR
- run the DEXSeq-based switch test when the design has two conditions with at least two samples each
- write notes, RDS objects, and switch tables when available

Reason for local implementation:

- this is the project-specific biological analysis step
- there is no generic nf-core module that exactly wraps this IsoformSwitchAnalyzeR use case

## Internal Data Flow

The downstream read-processing path uses tuples shaped like:

`tuple val(meta), path(reads)`

Important metadata fields include:

- `id`
- `condition`
- `replicate`
- `single_end`
- `strandedness`
- `batch`
- `sra_run`, the internal copy of the SRA manifest `run_accession` value

The original input metadata file is also passed to the ISAR module so sample-to-condition information remains available at cohort level.

## Implemented Milestones

- General FASTQ samplesheet parsing
- Multiple runs per sample via `cat/fastq`
- Minimal SRA manifest mode
- FastQC, fastp, Salmon, and MultiQC wiring
- Local IsoformSwitchAnalyzeR module
- Synthetic test profile with a tiny reference and 2-vs-2 design
- Real SRA validation on a small `GSE95132` 2-vs-2 subset
- Project metadata cleanup from template `anton/isoform` to `Anton-Bch/isoform-nf-core`

## Remaining Architecture Work

- Decide whether to add separate automated profiles for single-end FASTQ and SRA smoke tests.
- Decide whether to support more than one contrast in the ISAR wrapper.
- Decide whether `batch` should remain metadata-only or become part of an ISAR/statistical design.
- Add optional functional annotation modules only after the core reference-based path remains stable.

## Out Of Scope For V1

- De novo novel isoform discovery
- Dataset-branded workflow modes
- Hardcoded disease-specific condition names
- Mandatory Pfam/SignalP/IUPred2A/DeepLoc2/DeepTMHMM integration
