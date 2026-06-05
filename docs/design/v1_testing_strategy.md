# V1 Testing Status And Strategy

Last updated: 2026-06-02

This document describes the current validation state and the remaining testing work for V1.

## Current Automated Smoke Test

The `test` profile now uses tiny synthetic data stored in the repository:

- synthetic transcript FASTA
- synthetic genome FASTA
- synthetic GTF
- four tiny single-end FASTQ files
- 2 control and 2 treated samples
- a contrast file

The smoke test exercises:

- samplesheet parsing
- FastQC
- `cat/fastq`
- fastp
- Salmon index and quantification
- IsoformSwitchAnalyzeR import and switch testing
- MultiQC

The profile has been run successfully on the VM with Docker.

## Manual Real-Data Validation

Real-data validation is intentionally outside normal CI.

Validated so far:

- `GSE95132` minimal 1-vs-1 SRA run
  - confirmed SRA download/conversion, preprocessing, Salmon, ISAR import, and MultiQC
  - ISAR correctly skipped statistical testing because there was only one replicate per condition
- `GSE95132` 2-vs-2 SRA run
  - confirmed SRA mode, full downstream processing, and DEXSeq-based ISAR test execution
  - no genes passed the current switch cutoffs in that small subset, which is acceptable for a workflow validation run

These runs also revealed and fixed two useful real-world issues:

- SRA mode should use `prefetch` before `fasterq-dump`.
- The newest IsoformSwitchAnalyzeR BioContainer was missing `MASS`, so the ISAR module was pinned to a working BioContainer and the conda environment includes `r-mass`.

## What Is Still Missing

Recommended next tests:

- Add a dedicated paired-end FASTQ smoke test.
- Add a dedicated single-end FASTQ smoke test if the current test is kept focused on one layout.
- Add an SRA-mode automated test only if a tiny deterministic public accession or a stubbed approach is acceptable.
- Update `nf-test` snapshots once the output contract stabilizes.
- Add assertions for important ISAR files such as `analysis_notes.txt`, `switchAnalyzeRlist_imported.rds`, and switch tables when switches are detected.

## Testing Principles

- Use synthetic data for fast CI and deterministic snapshots.
- Use real public datasets for manual validation of runtime, SRA behavior, mapping rates, and biological plausibility.
- Avoid placing large public SRA downloads in routine CI.
- Prefer Docker or Singularity/Apptainer profiles for reproducibility.

## Current Minimum Acceptance

Before treating a future change as safe, run at least:

```bash
nextflow run . -profile test,docker --outdir results/test
```

For interface or SRA-related changes, also run a small real-data SRA validation on the VM.
