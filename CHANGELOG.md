# Anton-Bch/isoform-nf-core: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
