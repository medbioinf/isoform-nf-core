# V1 Architecture Plan

This document turns the V1 interface design into a concrete implementation plan for the `isoform-nf-core` pipeline.

It focuses on:

- which shared `nf-core` modules to use
- which local modules and subworkflows to create
- how the workflow should be split
- what should be implemented first

## Architectural goal

The future production pipeline should implement one clean, generic workflow:

`FASTQ or optional SRA manifest -> QC / trimming -> Salmon -> IsoformSwitchAnalyzeR -> MultiQC / reports`

The main design principle is:

- use shared `nf-core` components for common bioinformatics steps
- keep only genuinely project-specific logic local

## Top-level workflow design

The top-level workflow should have one clear branching point at the input layer:

1. FASTQ samplesheet input
2. SRA manifest input

After that branching point, both modes should converge into one shared internal sample representation and follow the same downstream path.

## Proposed workflow graph

### Input layer

Two mutually exclusive entry paths:

- `INPUT_FASTQ`
- `INPUT_SRA`

Both should emit:

- normalized `meta`
- normalized read files
- normalized cohort metadata file for downstream analysis

### Preprocessing layer

- `FASTQC`
- optional `FASTP`

This layer should take normalized reads and emit:

- QC reports
- cleaned reads

### Quantification layer

- `SALMON_INDEX`
- `SALMON_QUANT`

This layer should emit:

- Salmon index
- per-sample Salmon quantification directories

### Analysis layer

- `ISAR_ANALYSIS`

This layer should emit:

- ISAR result directory
- cohort-level result tables
- plots

### Reporting layer

- `MULTIQC`
- software versions collation
- workflow summary

## Shared nf-core modules to use

These modules should be reused from `nf-core` rather than reimplemented locally.

### Already present in the scaffold

- `fastqc`
- `multiqc`

### Recommended to install next

- `fastp`
- `salmon/index`
- `salmon/quant`

### Candidate to evaluate, but not mandatory for v1

- `sratools/prefetch`
- `fastq-dl`

These downloader-related components should only be adopted if they fit the actual SRA input design cleanly.

If they do not fit well, local SRA input modules are acceptable.

## Local modules to create

These parts are specific enough that they should likely remain local modules.

### `ISAR_ANALYSIS`

Purpose:

- take cohort metadata, quantification outputs, transcript FASTA, and GTF
- run `IsoformSwitchAnalyzeR`
- write result tables and plots

Why local:

- this is the project’s central biology-specific step
- it is not a standard shared `nf-core` module

### `WRITE_RUNTIME_SAMPLESHEET` or equivalent

Purpose:

- convert normalized in-memory metadata into a concrete CSV file for downstream tools when needed

Why local:

- this is workflow glue specific to your input model

### `FETCH_SRA_FASTQ` or split SRA modules

Purpose:

- download public sequencing data
- convert it to FASTQ
- normalize output naming

Why local:

- this behavior is strongly tied to your optional SRA input mode
- a shared downloader may not match your exact needs

Possible split:

- `SRA_PREFETCH`
- `SRA_FASTERQ_DUMP`
- optional compression / rename step

That split may be preferable if it improves retries and clarity.

## Local subworkflows to create

The following subworkflows should likely exist under `subworkflows/local/`.

### `input_fastq`

Responsibilities:

- validate FASTQ samplesheet rows
- normalize layout
- create the internal sample representation
- emit metadata CSV / channel for downstream use

### `input_sra`

Responsibilities:

- validate SRA manifest rows
- download / convert SRA data
- normalize output reads
- emit the same internal sample representation as FASTQ mode

### `preprocess_reads`

Responsibilities:

- run `FASTQC`
- optionally run `FASTP`
- return the reads that should go into Salmon

This should be written so the rest of the workflow does not care whether trimming was skipped.

### `salmon_quantification`

Responsibilities:

- build the decoy-aware Salmon index
- quantify each sample
- emit per-sample quantification outputs

### `isar_analysis`

Responsibilities:

- collect cohort-level quantification directories
- pass metadata and references to the local ISAR module
- emit result directory and key outputs

## Internal channel contract

The shared downstream workflow should operate on one normalized structure.

Suggested read channel shape:

- `tuple val(meta), path(reads)`

Suggested `meta` fields:

- `id`
- `condition`
- `layout`
- `patient_id`
- `sra_run`

Suggested metadata CSV shape:

- one normalized cohort-level CSV with the minimal columns required by the ISAR step

The important design rule is:

- FASTQ mode and SRA mode must converge into this same internal contract before preprocessing begins

## Recommended implementation sequence

### Milestone 1: schema and top-level contract

- update `nextflow.config`
- update `nextflow_schema.json`
- update `assets/schema_input.json` or replace it with the new samplesheet schema
- add schema support for `sra_manifest`

### Milestone 2: shared modules

- install `fastp`
- install `salmon/index`
- install `salmon/quant`

At this point the pipeline will have the core shared building blocks it needs.

### Milestone 3: local input subworkflows

- implement `input_fastq`
- implement `input_sra`
- normalize both to the same internal sample representation

### Milestone 4: local preprocessing and quantification wiring

- implement `preprocess_reads`
- implement `salmon_quantification`
- connect them to the top-level workflow

### Milestone 5: local ISAR module

- port the current reusable ISAR R script into the new repo
- wrap it as a local module with proper `nf-core` conventions
- wire it into `isar_analysis`

### Milestone 6: reporting and polishing

- connect `MultiQC`
- collate software versions
- update docs and citations
- add tests and run lint

## What should stay out of the first architecture

The following should not shape the first implementation:

- dataset-branded workflow branches
- CRC-specific validation gene panels as required defaults
- broad support for every optional ISAR companion tool
- many optional covariates not yet used in the workflow

These can be reconsidered later once the core pipeline is stable.

## Decision points to keep in mind

### Downloader implementation

Still open:

- use local `prefetch + fasterq-dump`
- or adopt a shared downloader approach if it fits

For v1, correctness and simplicity matter more than maximal reuse.

### ISAR packaging

Preferred first attempt:

- local module using package-based environment definitions

Fallback:

- one dedicated pinned custom container only for the ISAR module

### `patient_id`

For v1, preserve it in metadata and internal channels.

Do not let advanced matched-design modeling block the first implementation unless it is clearly required for correctness.

## Bottom line

The V1 architecture should be small and deliberate:

- shared `nf-core` modules for the standard steps
- local subworkflows for input normalization and project-specific orchestration
- one local ISAR module for the central analysis step

That gives the pipeline a strong `nf-core` backbone without forcing the unique isoform-switch logic into an awkward generic shape.
