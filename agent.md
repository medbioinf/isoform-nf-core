# Agent Guide

Last updated: 2026-04-22

## Purpose of this repo

This repository is the clean implementation target for the future `nf-core` isoform-switch pipeline.

Its job is to become a small, generic, reproducible pipeline for:

`FASTQ or optional SRA manifest -> QC / trimming -> Salmon -> IsoformSwitchAnalyzeR`

This repo should follow `nf-core` standards from the start.

The prototype and research workspace lives in the sibling repository:

- `nextflow-studienprojekt`

## Current status

This repo was created from the `nf-core` template and is still early.

Right now it contains:

- the standard `nf-core` scaffold
- initial planning docs for the V1 interface and architecture
- template modules already present for `fastqc` and `multiqc`

It does **not** yet contain the real workflow implementation.

## Key planning docs

Before making architecture changes or porting prototype logic, read:

- [docs/v1_interface_design.md](./docs/v1_interface_design.md)
  - current proposed V1 user-facing contract
- [docs/v1_architecture_plan.md](./docs/v1_architecture_plan.md)
  - current proposed module and subworkflow structure

These docs are the current design baseline.

## Core V1 scope

The intended V1 scope is:

- paired-end RNA-seq
- single-end RNA-seq
- primary input via FASTQ samplesheet
- optional input via SRA manifest
- transcript-level quantification with Salmon
- cohort-level isoform-switch analysis with `IsoformSwitchAnalyzeR`

V1 should not include:

- dataset-branded modes
- hardcoded GEO-specific defaults
- every optional external ISAR companion tool
- broad condition- or project-specific branching

## Current intended interface

### FASTQ samplesheet

Required columns:

- `sample`
- `condition`
- `fastq_1`

Optional columns:

- `fastq_2`
- `layout`
- `patient_id`

### SRA manifest

Required columns:

- `sample`
- `condition`
- `sra_run`

Optional columns:

- `layout`
- `patient_id`

### Required references

- `transcript_fasta`
- `genome_fasta`
- `gtf`

Do not reintroduce dataset-branded parameter modes unless there is a very strong reason.

## Architecture direction

The intended workflow structure is:

- input branch for FASTQ mode
- input branch for SRA mode
- shared preprocessing path
- shared Salmon quantification path
- local ISAR analysis path
- shared reporting path

Expected shared `nf-core` modules:

- `fastqc`
- `fastp`
- `multiqc`
- `salmon/index`
- `salmon/quant`

Expected local modules / subworkflows:

- input normalization
- optional SRA download / conversion logic
- runtime samplesheet writing if needed
- `IsoformSwitchAnalyzeR` analysis

## What to implement first

The recommended order is:

1. finalize `nextflow.config` params
2. update `nextflow_schema.json`
3. define the two input modes in code
4. install and wire shared `nf-core` modules
5. create the local ISAR module
6. add tests
7. run lint and iterate until clean

Do not start by copying the prototype repo wholesale.

## Standards that matter here

This repo should be developed to satisfy `nf-core` expectations:

- preserve template structure
- keep schema and params aligned
- use `nf-test`
- use pinned software environments
- prefer shared `nf-core` modules where possible
- keep local modules compliant with `nf-core` conventions
- keep the final pipeline generic, not dataset-branded

The detailed project-specific interpretation of these standards is documented in the prototype repo:

- `nextflow-studienprojekt/docs/nfcore_standards_for_this_project.md`

## Packaging strategy

Use shared module containers for standard tools.

Avoid one giant custom pipeline-wide Docker image.

For the ISAR step:

- first try a local module with package-based environment definitions
- only fall back to a dedicated custom container if package-based execution proves unstable

Do not use floating tags like `latest` in the final implementation.

## Testing strategy

This repo should use small synthetic data for CI and `nf-test`.

Minimum intended coverage:

- one paired-end smoke test
- one single-end smoke test
- stable minimal outputs

Real datasets such as `GSE50760` and `GSE95132` should be used for validation outside routine CI.

## Relationship to the prototype repo

Use `nextflow-studienprojekt` to learn:

- what metadata is truly needed
- what breaks on real cohorts
- which steps are stable enough to extract

Use this repo to implement:

- the clean general interface
- `nf-core` schema and docs
- the final workflow graph
- tests and lint-compliant structure

## Quick rule of thumb

If something is a standard bioinformatics building block, prefer shared `nf-core` modules.

If something is the unique biology-specific core of this project, keep it local but implement it in an `nf-core`-compatible way.
