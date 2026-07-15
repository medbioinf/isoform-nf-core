# Anton-Bch/isoform-nf-core: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0-beta.2 - 2026-07-15

Second beta release of Anton-Bch/isoform-nf-core, focused on scientific correctness, reproducibility, and release validation.

### `Added`

- Deterministic paired-end and multi-lane FASTQ fixtures with regression coverage for paired-read ordering.
- Multi-run SRA and per-isoform GO significance regression tests.
- Requested and inferred Salmon library-type provenance in pipeline outputs.

### `Changed`

- Differential isoform usage models now include batch and all complete additional covariates, with explicit confounding checks.
- Salmon automatic library inference remains automatic instead of being converted to an unstranded assumption.
- Optional annotation branches retain execution reports while tolerating containers without process-metric support.

### `Fixed`

- Preserved biological samples with multiple SRA runs through download, regrouping, concatenation, and paired-read ordering.
- Required GO q-value and dIF thresholds to pass on the same isoform, preventing cross-isoform false positives.
- Made Pfam indexing compatible with read-only reference directories and tolerated empty IUPred2A imports.
- Pinned DeepTMHMM and DeepLoc2 containers by immutable digest.
- Restored parser compatibility with Nextflow 26 by removing the deprecated `for` loop syntax.
- Hardened the privileged template-version workflow so pull-request code is never checked out or executed.
- Enabled release-mode nf-core linting for pull requests targeting `main` as well as `master`.

### `Dependencies`

- Requires Nextflow `>=25.04.0`.
- Uses nf-core template `3.5.2` conventions.
- Updates CI to nf-test `0.9.5` for Nextflow 26-compatible test wrappers.

## v1.0.0-beta.1 - 2026-07-07

First beta release of Anton-Bch/isoform-nf-core.

### `Added`

- Reference-based RNA-seq workflow from FASTQ or SRA manifest input through FastQC, fastp, Salmon, IsoformSwitchAnalyzeR, and MultiQC.
- Multi-contrast isoform switch summaries and lightweight IsoformSwitchAnalyzeR visualizations.
- Optional downstream GO enrichment and functional annotation support for Pfam, IUPred2A, SignalP, DeepTMHMM, DeepLoc2, and annotated switch plots.
- Synthetic test fixtures, nf-test coverage, and developer documentation for local validation.

### `Fixed`

- Hardened sample-sheet validation, samplesheet-relative FASTQ handling, and duplicate sample ID detection.
- Reduced CI resource profiles for hosted runner compatibility.
- Decoupled optional GO enrichment from the default ISAR runtime path.
- Aligned the full test profile with the current explicit reference-file contract and hosted-runner resource limits.

### `Dependencies`

- Requires Nextflow `>=25.04.0`.
- Uses nf-core template `3.5.2` conventions and containerized BioContainers-based workflow steps.

### `Deprecated`

- None.
