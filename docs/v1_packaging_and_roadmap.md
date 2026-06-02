# V1 Packaging And Roadmap Status

Last updated: 2026-06-02

This document summarizes the current packaging choices and the remaining roadmap after the first reusable workflow path was implemented.

## Packaging Choices

The pipeline currently follows this split:

- standard RNA-seq tools use nf-core modules and their normal Bioconda/BioContainer definitions
- project-specific logic lives in local modules

Standard module-backed tools:

- FastQC
- cat/fastq
- fastp
- Salmon
- MultiQC

Local project-specific modules:

- `FETCH_SRA_FASTQ`
- `ISAR_ANALYSIS`

## ISAR Runtime

The ISAR module is currently pinned to:

`quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0`

Reason:

- the newer `2.10.0` BioContainer was missing the R package `MASS`
- the `2.6.0` BioContainer successfully ran the real `GSE95132` 2-vs-2 validation
- the conda environment also includes `r-mass`

Future work:

- periodically re-check newer BioContainer builds
- consider a tiny project-owned ISAR container only if package-based environments remain fragile
- keep any custom container scoped to the ISAR module, not the whole pipeline

## SRA Runtime

SRA mode uses a local module with:

- `prefetch`
- `fasterq-dump`
- `pigz`

This was validated on real `GSE95132` runs. Direct `fasterq-dump` from accession was not robust enough in the container, so the module downloads the archive first and then converts the local file.

## Current Roadmap Status

Completed:

- Replace template identity and broken `anton/isoform` links.
- Implement FASTQ samplesheet input.
- Implement minimal SRA manifest input.
- Reuse nf-core modules for FastQC, cat/fastq, fastp, Salmon, and MultiQC.
- Add a local IsoformSwitchAnalyzeR module.
- Add synthetic test fixtures and a working `test,docker` profile.
- Validate the workflow on synthetic data.
- Validate real SRA execution on a small `GSE95132` 1-vs-1 and 2-vs-2 subset.

Still recommended before a polished release:

- Add stronger `nf-test` assertions and snapshots.
- Add separate single-end and paired-end test coverage if needed.
- Run `nf-core pipelines lint` and address template/lint findings.
- Improve the README with a workflow diagram.
- Expand `CITATIONS.md` for Salmon, fastp, IsoformSwitchAnalyzeR, DEXSeq, and SRA tools.
- Decide whether optional functional annotation modules belong in V1 or a later milestone.

## Optional Future Extensions

Not part of the current V1 path:

- de novo novel isoform discovery
- Pfam domain annotation
- SignalP
- IUPred2A
- DeepLoc2
- DeepTMHMM
- biological interpretation reports beyond the base ISAR outputs

These tools should be integrated only after the core reference-based workflow remains stable and the container/licensing situation is clear.

## Bottom Line

The project has moved beyond the original scaffold/planning phase. The current priority is no longer interface design, but hardening:

- tests
- linting
- citations
- docs polish
- optional annotation-module scoping
