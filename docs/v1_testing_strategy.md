# V1 Testing Strategy

This document defines the testing strategy for the first real implementation of `isoform-nf-core`.

The main goal is to make the pipeline:

- easy to validate during development
- compatible with `nf-core` expectations
- testable in CI without depending on large public datasets

## Testing principles

The pipeline should use two different kinds of validation:

### 1. CI and nf-test validation

Use small synthetic data and deterministic outputs.

Purpose:

- fast feedback
- stable `nf-test` snapshots
- compatibility with `nf-core` CI and lint expectations

### 2. Real-data validation

Use the existing prototype datasets outside CI.

Purpose:

- confirm behavior on real sequencing data
- check runtime, memory, and disk usage
- validate input assumptions and biological outputs

The important rule is:

- synthetic data is for automated tests
- real datasets are for manual or scheduled validation

## What should be tested in v1

The pipeline needs to support both read layouts and both input modes.

That means the testing plan should explicitly cover:

- paired-end FASTQ mode
- single-end FASTQ mode
- optional SRA mode
- the shared downstream path into Salmon
- the optional ISAR analysis step

## Recommended test layers

### Layer 1: Minimal pipeline smoke tests

These are the highest-priority tests.

They should confirm that the whole pipeline runs successfully on very small synthetic datasets.

Required smoke tests:

1. paired-end FASTQ smoke test
2. single-end FASTQ smoke test

These tests should:

- use tiny synthetic references
- use tiny synthetic FASTQ files
- complete quickly
- verify the basic workflow graph and file naming

### Layer 2: Optional ISAR-enabled smoke test

This test should run the downstream ISAR step on minimal synthetic data if runtime and determinism are acceptable.

Purpose:

- confirm the integration between Salmon output and the ISAR module
- catch breakage in the R environment or interface

This test does not need to prove biological realism.

It only needs to prove:

- the step runs
- outputs are created
- snapshots are stable enough for CI

### Layer 3: Optional SRA-mode smoke test

This should not be the first priority.

SRA mode adds download complexity and may be harder to keep deterministic in CI.

Recommended approach:

- first implement SRA mode
- validate it manually
- then decide whether a tiny mocked or lightweight automated test is practical

If a stable CI test is not practical at first, this can remain a manual validation path temporarily.

## Recommended V1 test matrix

### Required in early implementation

- `test`
  - paired-end FASTQ mode
  - `run_isar = false`
- `test_single`
  - single-end FASTQ mode
  - `run_isar = false`

These should be the first two test profiles brought into shape.

### Strongly recommended after the ISAR module is working

- `test_isar`
  - small synthetic data
  - `run_isar = true`

### Later / optional

- `test_sra`
  - optional SRA input path if reproducible enough for automation

## Synthetic test data requirements

The synthetic test data should be:

- very small
- deterministic
- stored specifically for this pipeline
- realistic enough to exercise both single-end and paired-end code paths

The data should include:

- tiny transcript FASTA
- tiny genome FASTA
- tiny GTF
- paired-end FASTQ example
- single-end FASTQ example
- matching samplesheets

The prototype repo already contains useful synthetic ideas that can be adapted, but the actual CI test assets should live in a stable form suitable for the `isoform-nf-core` pipeline.

## Snapshot philosophy

The tests should snapshot:

- stable output paths
- stable key text outputs
- software versions file, after removing the variable Nextflow version if needed

The tests should avoid snapshotting:

- unstable HTML output internals if they change too often
- timestamps
- run-specific trace metadata

This matches normal `nf-test` practice and keeps the tests maintainable.

## What to assert in the first tests

For the first pipeline smoke tests, useful assertions are:

- workflow succeeds
- expected output directories exist
- expected key files exist
- file names are stable

For example:

- `multiqc_report.html`
- Salmon quant directories
- expected `pipeline_info` outputs

For `test_isar`, additional useful assertions are:

- ISAR result directory exists
- key result tables exist

For example:

- `switch_summary.csv`
- `top_switches.csv`
- filtered annotation output if that remains part of the module design

## Real-data validation plan

Real-data validation should remain outside routine CI.

Use:

- `GSE50760` for paired-end real-data validation
- `GSE95132` for single-end and SRA-mode validation

These validation runs should answer:

- does the pipeline behave correctly on real cohorts?
- are the memory and runtime defaults reasonable?
- are there hidden assumptions that do not show up in synthetic tests?

## Immediate implications for the scaffold

The current scaffold still contains template placeholders that must be replaced.

Specifically:

- `conf/test.config`
  - still points at unrelated template test data
- `tests/nextflow.config`
  - still uses a placeholder test-data branch
- `tests/default.nf.test`
  - should be updated once the real V1 outputs are defined

These should be updated as soon as schema and input parsing are implemented.

## Recommended implementation order

1. build the synthetic paired-end FASTQ test data path
2. build the synthetic single-end FASTQ test data path
3. replace template `test` config with paired-end V1 test
4. add a `test_single` profile
5. update `tests/default.nf.test` to match the real outputs
6. add an ISAR-enabled minimal test once the local ISAR module is stable

## Bottom line

The V1 pipeline should be considered test-ready when it has:

- one paired-end smoke test
- one single-end smoke test
- deterministic synthetic references and reads
- `nf-test` snapshots for the stable output contract

That is the minimum testing baseline needed before the pipeline should be treated as a serious `nf-core` implementation candidate.
