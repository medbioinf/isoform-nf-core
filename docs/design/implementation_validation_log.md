# Implementation Validation Log

Last updated: 2026-06-12

This document records what we tested while implementing the current reusable pipeline. It complements the higher-level [V1 testing strategy](v1_testing_strategy.md) by documenting concrete validation runs, what each run was meant to prove, and what we learned.

The goal is not to claim that every step is fully production-tested. The goal is to keep a clear project record of how each implementation step was checked.

## Test Environments Used

### Local Development Machine

Used for:

- editing code and documentation
- checking git status and diffs
- copying selected result files back from the VM
- inspecting generated output tables and plots

Limitation:

- The local environment did not have Java available during the latest checks, so full Nextflow execution was validated on the VM instead.

### VM With Docker Runtime

Used for executable validation.

Typical profile:

```bash
-profile test,docker
```

Why the VM was used:

- Docker runtime was available.
- Bioinformatics containers could be pulled and executed.
- Real public RNA-seq data and reference resources were already staged there.
- Long-running steps such as SRA download and Pfam scan are more appropriate on the VM than on a laptop.

## Validation Principles

We used three levels of validation:

- Syntax/config validation: does Nextflow parse the workflow and schema?
- Synthetic smoke validation: does the workflow graph run on tiny test data?
- Real-data validation: does the implementation work with realistic public RNA-seq data and reference files?

Synthetic tests are fast and good for catching wiring mistakes. Real-data tests are slower but catch biological-data problems that tiny fixtures cannot represent.

## Core FASTQ Path

Implementation area:

- FASTQ samplesheet input
- read metadata parsing
- FastQC
- multiple-run concatenation with `cat/fastq`
- fastp trimming
- Salmon index and quantification
- MultiQC

Validation performed:

```bash
nextflow run . -profile test,docker --outdir results/test
```

What this checked:

- The test profile can run from repository fixtures without extra user data.
- FASTQ rows are parsed into the expected Nextflow `[meta, reads]` structure.
- Core nf-core modules execute inside containers.
- Salmon index and quantification complete on the synthetic reference.
- MultiQC collects reports and pipeline information.

Outcome:

- Passed on the VM.

Remaining gaps:

- Add dedicated small paired-end and single-end fixtures if we want explicit coverage for both layouts.
- Add nf-test assertions once the output contract stabilizes.

## SRA Input Path

Implementation area:

- minimal SRA manifest
- `prefetch`
- `fasterq-dump`
- gzip compression
- conversion into the same internal FASTQ channel as normal FASTQ input

Validation performed:

- Real-data SRA validation with a minimal `GSE95132` 1-vs-1 subset.
- Real-data SRA validation with a small `GSE95132` 2-vs-2 subset.

What this checked:

- The pipeline can start from public run accessions.
- SRA downloads are performed before conversion.
- Single-end versus paired-end output is detected after conversion.
- SRA-derived FASTQs continue through the same downstream path as local FASTQs.

Outcome:

- SRA download/conversion worked on the VM.
- The 1-vs-1 run correctly reached ISAR import but skipped statistical testing because each condition had only one replicate.
- The 2-vs-2 run reached DEXSeq-based ISAR testing.

Issues discovered and fixed:

- `fasterq-dump` alone was not robust enough for the intended SRA path; `prefetch` was added before conversion.

Remaining gaps:

- SRA mode still expects run accessions such as `SRR...`, `ERR...`, or `DRR...`.
- It does not yet resolve higher-level accessions such as `GSE...`, `GSM...`, or `SRP...`.
- Automated SRA testing is not yet included because it would depend on network access and public service availability.

## IsoformSwitchAnalyzeR Analysis

Implementation area:

- Salmon output import into IsoformSwitchAnalyzeR
- contrast handling
- DEXSeq-based switch test where enough replicates are available
- graceful handling of tiny smoke-test datasets

Validation performed:

- Synthetic `test,docker` run.
- Real-data `GSE95132` 1-vs-1 SRA run.
- Real-data `GSE95132` 2-vs-2 SRA run.
- Real-data `GSE50760` 4-vs-4 run used later for visualization and Pfam validation.

What this checked:

- ISAR can import Salmon quantification outputs.
- The GTF and transcript FASTA are accepted as reference inputs.
- The pipeline can distinguish import-only/smoke-test cases from cases where statistical switch testing is possible.
- ISAR output files are available for downstream visualization modules.

Outcome:

- ISAR import worked for synthetic and real-data runs.
- Statistical testing ran when sufficient replicates were available.
- The pipeline produced a reusable ISAR result directory for downstream visualization.

Issues discovered and fixed:

- The newest IsoformSwitchAnalyzeR BioContainer was missing the R package `MASS`.
- The ISAR module was pinned to a working BioContainer, and the conda fallback includes `r-mass`.

Remaining gaps:

- Add stronger automated assertions for key ISAR outputs.
- Add a stable fixture where at least one switch is expected, if feasible.

## ISAR Visualization

Implementation area:

- top switch plots
- summary plots and tables from ISAR output
- gene-level and isoform-level q-value explanations in visualization outputs
- robust handling when no significant switch is available

Validation performed:

- Synthetic `test,docker` run.
- Real-data `GSE50760` 4-vs-4 validation run.
- Copied selected visualization outputs back from the VM for local inspection.

What this checked:

- Visualization can run from the pipeline-produced ISAR directory.
- Plots are produced without requiring external annotation databases.
- Empty/no-switch cases exit successfully with notes instead of failing the workflow.
- Real-data top switch plots are generated from the nf-core pipeline output.

Outcome:

- The visualization module completed in the synthetic and real-data paths.
- The GSE50760 4-vs-4 run produced ISAR visualization outputs including top switch plots.

Issues discovered and fixed:

- The first page of the top-switch PDF was blank; this was fixed by adding a proper cover/intro page.
- A prototype-specific explanatory message was removed from the reusable pipeline output.

Remaining gaps:

- More visualization add-ons from the prototype still need to be transferred step by step.

## Pfam Annotation Path

Implementation area:

- extraction of amino-acid sequences from significant ISAR switch candidates
- running `pfam_scan.pl` inside the workflow when `--pfam_db` is provided
- importing Pfam hits into IsoformSwitchAnalyzeR with `analyzePFAM()`
- generating Pfam consequence summaries and domain architecture plots

Supported modes:

- Normal reproducible route: provide `--pfam_db`.
- Advanced/debug route: provide `--pfam_results` to import an existing `pfam_scan.pl` result.

### Pfam Config Check

Validation performed on the VM:

```bash
nextflow config -profile test,docker
```

What this checked:

- The workflow parses after adding Pfam modules and parameters.
- The schema/config additions do not break the test profile.

Outcome:

- Passed.

### Pfam Import Shortcut Check

Validation performed on the VM:

```bash
nextflow run . \
    -profile test,docker \
    -stub-run \
    --pfam_results results/gse50760_pfam_scan_smoke_validation/pfam_scan.out \
    --outdir results/test_pfam_import_after_default_tighten \
    -work-dir work_test_pfam_import_after_default_tighten
```

What this checked:

- The workflow can import an existing `pfam_scan.pl` result.
- `PFAM_IMPORT` and `PFAM_VISUALIZATION` are wired correctly.
- With the current defaults, `PFAM_PREPARE` is skipped when only `--pfam_results` is supplied.

Outcome:

- Passed.
- `PFAM_IMPORT`, `PFAM_VISUALIZATION`, and `MULTIQC` completed.

### Pfam Database Route Stub Check

Validation performed on the VM:

```bash
nextflow run . \
    -profile test,docker \
    -stub-run \
    --pfam_db ../nextflow-studienprojekt/gse50760/reference/pfam/Pfam37.0 \
    --outdir results/test_pfam_db_stub_full_route \
    -work-dir work_test_pfam_db_stub_full_route
```

What this checked:

- Providing `--pfam_db` automatically activates `PFAM_PREPARE`.
- `PFAM_PREPARE -> PFAM_SCAN -> PFAM_IMPORT -> PFAM_VISUALIZATION` is connected correctly in the Nextflow graph.
- The `PFAM_SCAN` process can be exercised in stub mode.

Outcome:

- Passed.
- `PFAM_PREPARE`, `PFAM_SCAN`, `PFAM_IMPORT`, `PFAM_VISUALIZATION`, and `MULTIQC` completed.

Issue discovered and fixed:

- The first `PFAM_PREPARE` stub emitted no FASTA, so `PFAM_SCAN` could not be exercised in `-stub-run`.
- The stub was updated to emit a tiny dummy amino-acid FASTA.

### Pfam Smoke Scan

Validation performed on the VM:

- Took three amino-acid sequences from a real GSE50760 Pfam candidate FASTA.
- Ran `pfam_scan.pl` in the Pfam BioContainer against the local Pfam database.
- Imported the resulting `pfam_scan.out` into ISAR.
- Ran Pfam visualization.

What this checked:

- The Pfam scan container can run against the local Pfam database.
- The generated `pfam_scan.pl` output format is accepted by `IsoformSwitchAnalyzeR::analyzePFAM()`.
- The downstream import and visualization scripts handle generated scan output.

Outcome:

- Passed.
- The smoke scan produced a valid `pfam_scan.out`.
- ISAR imported 9 domain rows for the tiny subset.
- No domain gain/loss consequences were found in the three-sequence subset, which is acceptable for a smoke test.

### Full GSE50760 Pfam Validation

Validation performed on the VM:

- Used the real GSE50760 4-vs-4 ISAR output produced by the reusable pipeline.
- Extracted amino-acid FASTA for significant switch candidates.
- Ran `pfam_scan.pl` against the local Pfam database.
- Imported the Pfam output into ISAR.
- Generated Pfam consequence and domain architecture visualizations.
- Copied selected non-RDS outputs back locally for inspection.

Key results:

- Significant switch candidates in ISAR feature table: `512`
- Top candidates written: `25`
- Amino-acid FASTA records written: `2919`
- Imported Pfam domain rows: `7019`
- Isoforms with `domain_identified == yes`: `2537`
- Domain consequence events: `134`
- Top domain architecture candidates plotted: `12`

Consequence summary:

```text
Domain loss, normal colon to primary crc: 39
Domain gain, primary crc to normal colon: 32
Domain gain, normal colon to primary crc: 31
Domain loss, primary crc to normal colon: 31
Domain switch, primary crc to normal colon: 1
```

What this checked:

- The Pfam route works with realistic ISAR output, not only synthetic data.
- The pipeline-produced ISAR object is compatible with sequence extraction and Pfam import.
- Domain consequence visualizations can be created from real Pfam annotations.
- Pfam scan is computationally non-trivial and should be treated as a real optional annotation module.

Outcome:

- Passed.

Performance note:

- The real Pfam scan/import path took noticeably longer than lightweight plotting steps.
- This supports keeping Pfam optional and requiring a prepared local Pfam database instead of downloading it automatically by default.

Copied-back local inspection path:

```text
results/gse50760_pfam_scan_direct_validation/
```

Large serialized R objects were intentionally not copied back from the VM.

## Documentation Validation

Implementation area:

- user-facing usage docs
- output docs
- interface explanations
- prerequisites document

Validation performed:

- Reviewed docs for consistency with implemented parameters.
- Updated Pfam docs after the interface changed from “external result import only” to “pipeline can run `pfam_scan.pl` with `--pfam_db`”.
- Added a dedicated prerequisites document.
- Linked the prerequisites page from the documentation index.

Outcome:

- Documentation now describes both the normal Pfam database route and the advanced Pfam result import route.
- The recommended user path is `--pfam_db`.
- `--pfam_results` is documented as advanced/debug reuse.

## Current Minimum Checks Before Future Commits

For normal implementation changes:

```bash
nextflow config -profile test,docker
nextflow run . -profile test,docker --outdir results/test
```

For changes touching Pfam workflow wiring:

```bash
nextflow run . \
    -profile test,docker \
    -stub-run \
    --pfam_db /path/to/prepared/pfam_db \
    --outdir results/test_pfam_db_stub
```

For changes touching Pfam import/visualization only:

```bash
nextflow run . \
    -profile test,docker \
    -stub-run \
    --pfam_results /path/to/pfam_scan.out \
    --outdir results/test_pfam_import_stub
```

For changes touching SRA mode:

- Run at least one small real SRA validation on the VM.
- Confirm `prefetch`, `fasterq-dump`, compression, and downstream FASTQ processing complete.

## Known Testing Limitations

- Full real-data runs are manual and not part of CI.
- SRA tests depend on network access and public archive availability.
- Pfam tests require a prepared local Pfam database.
- The synthetic fixture is intentionally tiny and cannot prove biological correctness.
- More nf-test snapshots/assertions should be added once the output contract stabilizes.
