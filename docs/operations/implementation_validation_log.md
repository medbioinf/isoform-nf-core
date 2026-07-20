# Implementation Validation Log

Last updated: 2026-07-20

This document records what we tested while implementing the current reusable pipeline. It documents concrete validation runs, what each run was meant to prove, and what we learned.

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

## IUPred2A Annotation Path

Implementation area:

- extraction of amino-acid sequences from significant ISAR switch candidates
- running IUPred2A with ANCHOR2 enabled
- converting multi-FASTA IUPred2A output into the block format imported by IsoformSwitchAnalyzeR
- importing IDR annotations with `analyzeIUPred2A()`

Implementation status:

- Added as an optional module controlled by `--run_iupred2a`.
- Added dedicated preparation, run, and import modules.
- Added a converter script based on the GSE50760 prototype.
- Added output, interface, and prerequisites documentation.

Validation performed locally:

- Python syntax check for `convert_iupred2a_to_isar.py`.
- Tiny converter smoke test with two artificial FASTA records.
- VM `nextflow config -profile test,docker`.
- VM `nextflow run . -profile test,docker -stub-run --run_iupred2a`.
- Direct VM container smoke test with a tiny amino-acid FASTA.

Real-data subset validation performed on the VM:

- Used 10 real GSE50760 isoforms selected from the Pfam candidate FASTA.
- Ran the IUPred2A container with ANCHOR2 enabled.
- Converted the raw multi-FASTA output into the IsoformSwitchAnalyzeR import format.
- Imported the converted output into the real GSE50760 ISAR object.

Outcome:

- Converter restored per-isoform output blocks correctly.
- Nextflow stub workflow completed successfully.
- The IUPred2A container produced IUPred2/ANCHOR2 rows for the tiny FASTA.
- The real-data subset import produced 34 IDR rows.
- IsoformSwitchAnalyzeR reported IDR information added to 8 transcripts.

Runtime note:

- The current third-party IUPred2A image lacks `ps`, which Nextflow needs for trace/timeline/report metrics.
- The pipeline now accepts `--iupred2a_container` so sites can provide a compatible image without changing module code.
- The NFkB real-data validation below used a VM-local derivative containing `procps`; a portable replacement image still needs to expose `python3`, `/opt/iupred2a/iupred2a.py`, and `ps`.

## SignalP Annotation Path

Implementation area:

- extraction of amino-acid sequences from significant ISAR switch candidates
- running SignalP 5 in eukaryotic mode
- importing signal peptide annotations with `analyzeSignalP()`

Implementation status:

- Added as an optional module controlled by `--run_signalp`.
- Added dedicated preparation, run, and import modules.
- Added output, interface, and prerequisites documentation.
- Replaced the unaffiliated public image dependency with a required user-provided, appropriately licensed SignalP 5.0b container.

Validation performed:

- VM `nextflow config -profile test,docker`.
- VM `nextflow run . -profile test,docker -stub-run --run_signalp`.
- Direct VM container smoke test with a tiny amino-acid FASTA.

Real-data subset validation performed on the VM:

- Used 10 real GSE50760 isoforms selected from the Pfam candidate FASTA.
- Ran the SignalP container after confirming the module-level `PATH` fix.
- Imported the SignalP output into the real GSE50760 ISAR object.

Full real-data workflow validation performed on the VM:

```bash
nextflow run . -profile docker \
    --input assets/gse50760_4v4_existing_fastp_samplesheet.csv \
    --contrasts assets/gse50760_primary_crc_vs_normal_colon.csv \
    --transcript_fasta ../nextflow-studienprojekt/gse50760/reference/gencode.v49.transcripts.fa.gz \
    --gtf ../nextflow-studienprojekt/gse50760/reference/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz \
    --run_signalp \
    --signalp_container <licensed-signalp-5.0b-image-or-sif> \
    --outdir results/gse50760_4v4_signalp_full_validation \
    --multiqc_title GSE50760_4v4_signalp_full_validation \
    -work-dir work_gse50760_4v4_signalp_full_validation
```

The validation samplesheet reused existing GSE50760 fastp FASTQs on the VM to avoid repeating SRA download and trimming. This still exercised the real Salmon, ISAR, SignalP preparation, SignalP run, and SignalP import path.

Outcome:

- Nextflow stub workflow completed successfully.
- The SignalP container produced short-format SignalP 5 output for the tiny FASTA.
- The real-data subset predicted and imported 1 signal peptide from 10 isoforms.
- The full real-data workflow completed successfully in 52m 9s with 39 successful tasks.
- SignalP preparation found 516 significant switch candidates and wrote 2924 amino-acid FASTA records.
- `signalp5_summary.signalp5` contained 2924 isoform predictions plus the SignalP header lines.
- SignalP import added 343 signal peptide rows.
- `isoform_features` contained 343 isoforms with `signal_peptide_identified == yes`, 69631 with `no`, and 35701 with `unknown`.

Issue discovered and fixed:

- Direct module execution initially failed with `signalp: command not found` because `bash -lc` reset the container `PATH`.
- Calling `/opt/signalp/bin/signalp` directly was not suitable because SignalP's bundled assets are resolved relative to its expected runtime layout.
- The module now exports `/opt/signalp/bin` onto `PATH` before invoking `signalp`.

Release note:

- The historical VM validation used a third-party SignalP image. That image is no longer pulled or redistributed by the release pipeline.
- Current runs must provide `--signalp_container` with a SignalP 5.0b image or SIF provisioned under the user's institutional license. Container contents and metric tooling are therefore user-controlled.

Copied-back local inspection path:

```text
results/gse50760_4v4_signalp_full_validation_copy/
```

The copied files include selected SignalP import outputs, raw SignalP output, pipeline parameters, and the DAG. They are local inspection artifacts and remain outside version control.

## DeepTMHMM, DeepLoc2, and Annotated Switch Plots

Implementation area:

- extraction of amino-acid sequences from significant ISAR switch candidates
- running DeepTMHMM and converting its three-line topology output for IsoformSwitchAnalyzeR
- running DeepLoc2 and converting its localization table for IsoformSwitchAnalyzeR
- chaining optional annotation imports so later imports preserve earlier annotation layers
- creating annotated `switchPlot()` outputs from the newest available annotated ISAR object

Implementation status:

- Added optional modules controlled by `--run_deeptmhmm`, `--run_deeploc2`, and `--run_annotated_switch_plots`.
- Added dedicated preparation, run, and import modules for DeepTMHMM and DeepLoc2.
- Added converter scripts based on the GSE50760 prototype output formats.
- Added annotated switch plots with automatic annotation-layer status reporting.
- Updated the annotation import scripts so Pfam, IUPred2A, SignalP, DeepTMHMM, and DeepLoc2 can be chained without discarding previous annotations.

Validation performed locally:

- Python syntax check for the DeepTMHMM and DeepLoc2 converter scripts.
- `git diff --check`.

Validation performed on the VM:

```bash
nextflow config -profile test,docker
nextflow run . -profile test,docker -stub-run -resume \
    --run_deeptmhmm true \
    --run_deeploc2 true \
    --deeploc2_container <licensed-deeploc-2.1-image-or-sif> \
    --run_annotated_switch_plots true \
    --outdir results/test_milestone_b_stub3 \
    -work-dir work_test_milestone_b_stub3
```

After the documentation refresh, the same graph was sanity-checked again on the VM with:

```bash
nextflow config -profile test,docker
nextflow run . -profile test,docker -stub-run -resume \
    --run_deeptmhmm \
    --run_deeploc2 \
    --deeploc2_container <licensed-deeploc-2.1-image-or-sif> \
    --run_annotated_switch_plots \
    --outdir results/test_milestone_b_stub_final \
    -work-dir work_test_milestone_b_stub_final
```

Outcome:

- Passed.
- `DEEPTMHMM_PREPARE`, `DEEPTMHMM_RUN`, `DEEPTMHMM_IMPORT`, `DEEPLOC2_PREPARE`, `DEEPLOC2_RUN`, `DEEPLOC2_IMPORT`, `ANNOTATED_SWITCH_PLOTS`, and `MULTIQC` completed in stub mode.

Full real-data validation performed on the VM:

```bash
nextflow run . -profile docker -resume \
    --input assets/gse50760_4v4_existing_fastp_samplesheet.csv \
    --contrasts assets/gse50760_primary_crc_vs_normal_colon.csv \
    --transcript_fasta ../nextflow-studienprojekt/gse50760/reference/gencode.v49.transcripts.fa.gz \
    --gtf ../nextflow-studienprojekt/gse50760/reference/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz \
    --run_deeptmhmm true \
    --run_deeploc2 true \
    --deeploc2_container <licensed-deeploc-2.1-image-or-sif> \
    --run_annotated_switch_plots true \
    --annotated_switch_top_n 5 \
    --outdir results/gse50760_4v4_milestone_b_full_validation \
    --multiqc_title GSE50760_4v4_milestone_b_full_validation \
    -work-dir work_gse50760_4v4_milestone_b_full_validation
```

The validation samplesheet reused existing GSE50760 fastp FASTQs on the VM to avoid repeating SRA download and trimming. This still exercised the real Salmon, ISAR, DeepTMHMM preparation/run/import, DeepLoc2 preparation/run/import, annotated switch plots, and MultiQC path.

Key results:

- Full run completed successfully on 2026-06-21.
- GSE50760 4-vs-4 ISAR analysis found 521 significant isoform candidates, 447 switches, and 476 switching genes for `normal_colon vs primary_crc`.
- DeepTMHMM preparation wrote 2958 amino-acid FASTA records.
- DeepTMHMM produced `predicted_topologies.3line` and `deeptmhmm_regions_isar.tsv`; the converted topology file contained 7049 imported topology rows.
- DeepLoc2 preparation wrote 2958 amino-acid FASTA records.
- DeepLoc2 produced `results_20260621-113336.csv` and `deeploc2_isar.csv`; DeepLoc2 import added `sub_cell_location` values for 2842 isoforms.
- Annotated switch plots were generated for ATP8B1, AMFR, EPB41L1, TEX2, and DSC2.
- `annotation_status.csv` reported ORF/PTC, DeepLoc2, and DeepTMHMM as available; Pfam, SignalP, and IUPred2A were correctly reported as missing because they were not enabled in this validation run.

Issues discovered and fixed:

- The DeepLoc2 image has an entrypoint that interfered with Nextflow command execution and even stub commands. The module now clears the container entrypoint and sets writable `HOME` and `MPLCONFIGDIR` paths.
- DeepTMHMM could not run directly in `/openprotein` because it needs writable working directories. The module now creates a writable work directory and symlinks the model/runtime assets.
- DeepTMHMM's bundled `predict.py` wrote one report to `/deeptmhmm_results.md`, which is not writable in the container. The module patches that path to a relative filename inside the writable work directory before execution.
- The first DeepTMHMM workdir symlink attempt excluded package directories and triggered `No module named 'experiments'`; the module now symlinks all required `/openprotein` entries except `predict.py`, which is copied and patched.
- A later prototype-parity check showed that explicit annotated switch plot genes can be absent from the automatic annotation target set when upstream switch statistics differ. The prepare modules now force-include `--annotated_switch_genes` in enabled annotation FASTA preparation without changing the original ISAR switch statistics.

Runtime notes:

- In this validation run, DeepLoc2 processed 2958 proteins in about 33 minutes and downloaded ESM/model assets on first use.
- DeepTMHMM processed 2958 proteins with a long embedding phase followed by topology prediction; it dominated the optional annotation runtime.
- Because these tools are heavy, they should remain optional interpretation modules.

Copied-back local inspection path:

```text
results/gse50760_4v4_milestone_b_full_validation/
```

The copied files include selected DeepTMHMM, DeepLoc2, annotated switch plot, ISAR, MultiQC, and pipeline parameter outputs. Large serialized R objects and heavy work-directory intermediates were intentionally left on the VM.

## Target Gene Annotation Preparation Validation

Validation command, run on the VM:

```bash
nextflow run . \
    -profile docker \
    -resume \
    --input assets/gse50760_4v4_existing_fastp_samplesheet.csv \
    --contrasts assets/gse50760_primary_crc_vs_normal_colon.csv \
    --transcript_fasta ../nextflow-studienprojekt/gse50760/reference/gencode.v49.transcripts.fa.gz \
    --gtf ../nextflow-studienprojekt/gse50760/reference/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz \
    --run_pfam_prepare \
    --annotated_switch_genes ZNRF3,PBX3,YEATS4 \
    --outdir results/gse50760_4v4_target_gene_prepare_validation5 \
    --multiqc_title GSE50760_4v4_target_gene_prepare_validation5 \
    -work-dir work_gse50760_4v4_target_gene_prepare_validation5
```

Outcome:

- The run completed successfully in 4m 13s with 99% cached work; only `PFAM_PREPARE` and `MULTIQC` ran fresh.
- `pfam_prepare_notes.txt` reported `Target genes requested: ZNRF3, PBX3, YEATS4`, `Isoform rows force-included by target gene: 17`, and `AA FASTA records written: 2968`.
- The target-gene candidate table contains `selection_reason=target_gene` rows for ZNRF3, PBX3, and YEATS4.
- The AA FASTA contains target protein-coding isoforms including `ENST00000406323.3` for ZNRF3, `ENST00000342287.9` for PBX3, and `ENST00000247843.7`, `ENST00000548020.5`, `ENST00000549685.5`, and `ENST00000552955.1` for YEATS4. The retained-intron YEATS4 transcript does not emit an AA sequence, as expected.

Issues discovered and fixed:

- After normalizing CLI option names from hyphens to underscores, the prepare scripts still checked for `isar-dir`. They now check and read `isar_dir`.
- Target isoforms below the dIF cutoff were initially set to exactly the cutoff for sequence extraction. ISAR's switching-gene extraction behaves as a strict cutoff, so genes such as YEATS4 with all target isoforms below the cutoff were still excluded from AA FASTA output. The prepare scripts now set target-only extraction dIF just above the cutoff while preserving the original reported statistics in the CSV and ISAR result.

Copied-back local inspection path:

```text
results/gse50760_4v4_target_gene_prepare_validation5/
```

## Target Gene Pfam Full Validation

Validation command, run on the VM:

```bash
nextflow run . \
    -profile docker \
    -resume \
    --input assets/gse50760_4v4_existing_fastp_samplesheet.csv \
    --contrasts assets/gse50760_primary_crc_vs_normal_colon.csv \
    --transcript_fasta ../nextflow-studienprojekt/gse50760/reference/gencode.v49.transcripts.fa.gz \
    --gtf ../nextflow-studienprojekt/gse50760/reference/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz \
    --pfam_db ../nextflow-studienprojekt/gse50760/reference/pfam/Pfam37.0 \
    --run_annotated_switch_plots \
    --annotated_switch_genes ZNRF3,PBX3,YEATS4 \
    --annotated_switch_condition1 normal_colon \
    --annotated_switch_condition2 primary_crc \
    --outdir results/gse50760_4v4_target_gene_pfam_full_validation \
    --multiqc_title GSE50760_4v4_target_gene_pfam_full_validation \
    -work-dir work_gse50760_4v4_target_gene_pfam_full_validation
```

Outcome:

- The run completed successfully in 28m 7s with 81.6% cached work.
- Fresh processes: `PFAM_PREPARE`, `PFAM_SCAN`, `PFAM_IMPORT`, `PFAM_VISUALIZATION`, `ANNOTATED_SWITCH_PLOTS`, and `MULTIQC`.
- `PFAM_SCAN` used the local indexed Pfam37.0 database and scanned the target-gene-expanded FASTA.
- `PFAM_IMPORT` imported 7211 domain rows and reported 2614 isoforms with `domain_identified == yes`.
- Annotated switch plots were generated for ZNRF3, PBX3, and YEATS4.
- `annotation_status.csv` reports ORF/PTC and Pfam protein domains as available. SignalP, IUPred2A, DeepLoc2, and DeepTMHMM are missing in this run because only Pfam was enabled.

Target gene checks:

- ZNRF3 plot contains `ZNRF_3_ecto` and `zf-RING_2` domain tracks across all three plotted isoforms.
- PBX3 plot contains `PBC` and `Homeobox_KN` domain tracks.
- YEATS4 plot contains `YEATS` domain tracks for the protein-coding isoforms.
- The retained-intron YEATS4 transcript remains unannotated for Pfam, as expected.

Copied-back local inspection path:

```text
results/gse50760_4v4_target_gene_pfam_full_validation/
```

## Multi-Contrast ISAR Summary Validation

Implementation area:

- contrast-aware ISAR statistical-test guard
- cross-contrast significant-switch summary
- UpSet-style isoform switch intersection plot

Validation performed:

- Ran `nextflow config -profile test,docker` on the VM.
- Ran a VM stub workflow:

```bash
nextflow run . \
    -profile test,docker \
    -stub-run \
    --outdir results/test_isar_contrast_summary_stub \
    -work-dir work_test_isar_contrast_summary_stub
```

- Ran a real tiny VM workflow:

```bash
nextflow run . \
    -profile test,docker \
    --outdir results/test_isar_contrast_summary_real \
    -work-dir work_test_isar_contrast_summary_real
```

- Ran `run_isar_contrast_summary.R` directly against a synthetic analyzed ISAR-like object with three comparisons to exercise multi-set intersection plotting.

Outcome:

- The stub workflow completed and instantiated `ISAR_CONTRAST_SUMMARY`.
- The real tiny workflow completed and wrote `isar/isar_analysis/comparisons.csv` plus `isar/isar_contrast_summary/`.
- The tiny workflow had one comparison, so it wrote the per-comparison bar plot and skipped the UpSet plot with an explanatory note.
- The synthetic multi-comparison check wrote `significant_isoform_switches_per_comparison.*`, `isoform_switch_intersections.csv`, `isoform_switch_intersection_members.csv`, and `isoform_switch_upset.*`.
- The synthetic intersections matched the expected memberships: one B-only intersection, one A+B overlap, one A+C overlap, and one C-only intersection.

Known limitation:

- A direct standalone attempt to fake a four-condition ISAR analysis from duplicated tiny fixture quantifications failed inside `IsoformSwitchAnalyzeR::importRdata()` before statistical testing. This was treated as a fixture limitation rather than a pipeline regression; the real pipeline test and direct summary-script test passed.

## Full Annotation Parity Validation

Validation command, run detached on the VM:

```bash
nextflow run . \
    -profile docker \
    -resume \
    --input assets/gse50760_4v4_existing_fastp_samplesheet.csv \
    --contrasts assets/gse50760_primary_crc_vs_normal_colon.csv \
    --transcript_fasta ../nextflow-studienprojekt/gse50760/reference/gencode.v49.transcripts.fa.gz \
    --gtf ../nextflow-studienprojekt/gse50760/reference/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz \
    --pfam_db ../nextflow-studienprojekt/gse50760/reference/pfam/Pfam37.0 \
    --run_deeptmhmm \
    --run_deeploc2 \
    --deeploc2_container <licensed-deeploc-2.1-image-or-sif> \
    --run_annotated_switch_plots \
    --annotated_switch_genes ZNRF3,PBX3,YEATS4 \
    --annotated_switch_condition1 normal_colon \
    --annotated_switch_condition2 primary_crc \
    --outdir results/gse50760_4v4_parity_overnight_20260622 \
    --multiqc_title GSE50760_4v4_parity_overnight_20260622 \
    -work-dir work_gse50760_4v4_parity_overnight_20260622
```

Outcome:

- The run completed successfully in 2h 55m 45s with 48 succeeded processes.
- `annotation_status.csv` reports ORF/PTC, Pfam protein domains, DeepLoc2 subcellular locations, and DeepTMHMM topology as available. SignalP and IUPred2A are missing because they were not enabled in this run.
- Annotated switch plots were generated for ZNRF3, PBX3, and YEATS4.
- Visual inspection confirmed that the ZNRF3 plot contains Pfam domains, DeepLoc2 locations, and DeepTMHMM topology tracks.
- PBX3 and YEATS4 plots contain their expected Pfam domain tracks and localization labels; DeepTMHMM topology tracks are rendered where DeepTMHMM predicts regions.
- The annotated plots were copied back locally to `results/gse50760_4v4_parity_overnight_20260622_annotated/annotated_switch_plots/`.

Issue discovered and fixed:

- `deeptmhmm_import_notes.txt` incorrectly reported `Isoforms with topology_identified == yes: 0` because the import script expected a `topology_identified` column in `isoformFeatures`. In the tested IsoformSwitchAnalyzeR version, `analyzeDeepTMHMM()` stores topology rows in `topologyAnalysis` and `switchPlot()` renders from that table. The import script now derives `topology_identified` from unique `topologyAnalysis$isoform_id` values.
- Focused VM validation against the completed parity run wrote corrected DeepTMHMM import outputs: 7131 topology rows, 2948 unique isoforms with topology rows, and 2948 `topology_identified == yes` isoforms. Among the requested target genes, 13 of 17 isoform rows have topology annotations.

## Full All-Annotation Parity Validation

Validation command, run detached on the VM:

```bash
nextflow run . \
    -profile docker \
    --sra_manifest assets/gse50760_4v4_sra_manifest.csv \
    --contrasts assets/gse50760_primary_crc_vs_normal_colon.csv \
    --transcript_fasta ../nextflow-studienprojekt/gse50760/reference/gencode.v49.transcripts.fa.gz \
    --gtf ../nextflow-studienprojekt/gse50760/reference/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz \
    --pfam_db ../nextflow-studienprojekt/gse50760/reference/pfam/Pfam37.0 \
    --run_signalp \
    --signalp_container <licensed-signalp-5.0b-image-or-sif> \
    --run_iupred2a \
    --run_deeptmhmm \
    --run_deeploc2 \
    --deeploc2_container <licensed-deeploc-2.1-image-or-sif> \
    --run_annotated_switch_plots \
    --annotated_switch_genes ZNRF3,PBX3,YEATS4 \
    --annotated_switch_condition1 normal_colon \
    --annotated_switch_condition2 primary_crc \
    --outdir results/20260623_gse50760_4v4_all_annotations_parity \
    --multiqc_title 20260623_gse50760_4v4_all_annotations_parity \
    -work-dir work_20260623_gse50760_4v4_all_annotations_parity
```

Outcome:

- The run completed successfully in 3h 45m 25s with 62 succeeded processes.
- `annotation_status.csv` reports ORF/PTC, Pfam protein domains, SignalP signal peptide, IUPred2A/NetSurfP IDR, DeepLoc2 subcellular locations, and DeepTMHMM topology as available.
- Annotated switch plots were generated for ZNRF3, PBX3, and YEATS4 with all enabled annotation layers available to `switchPlot()`.
- SignalP imported 353 signal peptide rows and reported 353 `signal_peptide_identified == yes` isoforms.
- DeepTMHMM imported 7117 topology rows and reported 2997 unique isoforms with topology rows.
- IUPred2A imported 4044 IDR rows.
- The full VM result directory was about 142 GB because it includes SRA/FASTQ intermediates. A curated local inspection subset was copied to `results/20260623_gse50760_4v4_all_annotations_parity/` in this repository. The local subset is about 4.1 GB and includes annotated switch plots, ISAR outputs, annotation import outputs, summary plots, MultiQC, pipeline info, and relevant RDS files, but intentionally excludes large raw/intermediate `fetch/`, `cat/`, `fastp/`, and full Salmon output directories.

Issue discovered and fixed:

- `iupred2a_import_notes.txt` incorrectly reported `Isoforms with IDR_identified == yes: 0` because the import script expected an `IDR_identified` column in `isoformFeatures`. In the tested IsoformSwitchAnalyzeR version, `analyzeIUPred2A()` stores IDR rows in `idrAnalysis` and `switchPlot()` renders from that table. The import script now derives `IDR_identified` from unique `idrAnalysis$isoform_id` values.
- Focused VM validation against the completed all-annotation run wrote corrected IUPred2A import outputs: 4044 IDR rows, 1549 unique isoforms with IDR rows, and 1549 `IDR_identified == yes` isoforms. Among the requested target genes, 8 of 17 isoform rows have IDR annotations.

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

## GO Enrichment CSV Handoff Validation

Date: 2026-07-07

Implementation area:

- optional GO enrichment module
- `bin/run_go_enrichment.R`
- `bin/run_isar_analysis.R`
- `modules/local/isoform_go_enrichment`

Issue discovered:

- The first merged GO enrichment implementation read `switchAnalyzeRlist_analyzed.rds` directly inside the GO enrichment module.
- The configured GO container, `docker.io/jungwooseok/webgestalt:1.0.3`, contains `WebGestaltR` but not `IsoformSwitchAnalyzeR`.
- A direct VM validation against the completed GSE50760 ISAR result failed before enrichment with:

```text
unable to find required package 'IsoformSwitchAnalyzeR'
there is no package called 'IsoformSwitchAnalyzeR'
```

Design change:

- `run_isar_analysis.R` now exports `isar/isar_analysis/go_gene_scores.csv`.
- `run_go_enrichment.R` consumes `--gene-score-file` instead of reading the serialized ISAR R object.
- `ISOFORM_GO_ENRICHMENT` passes `isar_results/go_gene_scores.csv` to the GO script.
- This keeps the GO module gene-level and WebGestalt-focused. It also removes the `IsoformSwitchAnalyzeR` runtime dependency from the GO enrichment container.

VM validation performed:

```bash
nextflow run . \
    -profile docker,test \
    --run_go_enrichment \
    --outdir results/20260707_go_enrichment_csv_handoff_test \
    -work-dir work_20260707_go_enrichment_csv_handoff_test
```

Outcome:

- The small real Docker test completed successfully.
- `isar/isar_analysis/go_gene_scores.csv` was produced.
- `go/go_input/go_gene_scores.csv` was produced as the GO-module input copy.
- `go/go_enrichment.csv` and `go/go_summary.txt` were produced.
- The test dataset had only two background genes, so WebGestaltR was correctly skipped with a clear summary note rather than failing.

Additional real-data validation:

- Used the completed GSE50760 4-vs-4-vs-4 ISAR result from `results/20260705_gse50760_4v4v4_all_annotations/isar/isar_analysis/`.
- Exported a real GO gene-score table with 40,464 contrast-gene rows.
- Ran `bin/run_go_enrichment.R` in the WebGestalt-only container using the exported CSV.

Outcome:

- GO enrichment ran successfully without `IsoformSwitchAnalyzeR` in the GO container.
- Three contrasts were evaluated.
- Significant gene counts:
  - `liver_metastasis_vs_normal_colon`: 468
  - `primary_crc_vs_normal_colon`: 729
  - `liver_metastasis_vs_primary_crc`: 168
- Two enriched biological-process terms were reported for `liver_metastasis_vs_primary_crc`.
- The other two contrasts completed cleanly and reported no enriched terms at FDR 0.05.

Operational note:

- A broader GSE50760 `-resume` attempt was stopped because Nextflow began scheduling upstream SRA fetch tasks again. This was not needed to validate the GO runtime fix and would have wasted VM time.
- The validation status after this change is: GO enrichment is real-data script/runtime validated and full Nextflow-wiring validated on the test dataset. A raw-SRA-to-GO full GSE50760 run after this refactor remains optional, expensive, and not required for the container dependency fix.

## NFkB Mouse Study With Precomputed Salmon Input

Date: 2026-07-17, corrected plot rerun on 2026-07-20

Implementation areas:

- new `--salmon_input` mode for existing Salmon quantification directories
- four contrasts in one ISAR run
- mouse GO enrichment
- Pfam, IUPred2A, and DeepTMHMM annotation
- annotated switch plots and MultiQC

Dataset and setup:

- Mouse macrophages from wild type and three NFkB-inhibitor knockout genotypes: Bcl3 (`B3KO`), IkBNS (`NSKO`), and IkBz (`ZKO`).
- The complete study contained controls, 3-hour infection, and 24-hour infection samples. This validation used the 24 samples needed for four control-versus-3-hour comparisons, with three replicates per condition.
- Existing `quant.sf` directories were supplied with `--salmon_input`; the read-QC, trimming, indexing, and quantification stages were intentionally skipped.
- Matching Ensembl mouse release 114 transcript FASTA and GTF files were supplied.
- SignalP and DeepLoc2 were not enabled because compatible licensed containers were not available for this run.

VM run:

```text
Repository: /home/ubuntu/projects/isoform-nf-core-nfkb-validation
Results: results/20260717_nfkb_3h_precomputed_salmon_annotations
Work: work_20260717_nfkb_3h_precomputed_salmon_annotations
```

Outcome:

- Pipeline completed successfully in 2 h 51 min: 9 processes succeeded and 7 were restored from cache.
- ISAR evaluated all four requested contrasts and reported 3144 switches across 2529 genes in the combined summary.
- GO enrichment evaluated 31,321 contrast-gene rows. Five biological-process terms passed FDR 0.05 for `B3KO_control_vs_B3KO_3h`; the other contrasts completed successfully with no terms passing the cutoff.
- Pfam scanning and import completed, IUPred2A imported non-empty IDR annotations, and DeepTMHMM imported non-empty topology annotations.
- Annotated switch plots reported ORF/PTC, Pfam, IUPred2A, and DeepTMHMM layers as available.
- A curated 26 MB inspection package was copied locally to `results/20260717_nfkb_3h_precomputed_salmon_annotations_interpretable/`; large RDS objects and raw predictor outputs were intentionally omitted.

Issues discovered and fixed:

- `quant_dir` was initially treated as a possible statistical covariate. It is now explicitly classified as a technical input column.
- Version suffixes on Ensembl transcript identifiers could prevent matching between Salmon and the annotation. ISAR import now uses `ignoreAfterPeriod = TRUE` in addition to `ignoreAfterBar = TRUE`.
- WebGestalt identifier errors were difficult to diagnose. GO enrichment now validates the selected identifier type against the selected organism before submitting enrichment calls.
- The default IUPred2A image lacks `ps`. The new `--iupred2a_container` parameter allowed validation with a compatible VM-local image containing `procps`.
- Five annotated switch plots were initially emitted as `NA` placeholders although their gene mappings were present. R logical filters had allowed unrelated rows containing missing values into each gene subset. Comparison and gene matching are now explicitly missing-value-safe, and automatic ranking excludes candidates without usable identifiers.

Focused corrected-plot validation:

- Resumed the completed workflow twice; all expensive annotation processes remained cached.
- All ten annotated switch plots were generated with valid names and non-empty content.
- The corrected set includes `Tmem164`, `Syk`, `Zfx`, `Nfkb2`, and `Sumf1` instead of five `NA_annotated_switch` placeholders.
- Visual inspection confirmed transcript structures, expression/usage panels, Pfam or IUPred tracks where present, and DeepTMHMM topology tracks.

Precomputed-input wiring check on 2026-07-20:

- `nextflow config -profile test,docker` parsed successfully on the VM.
- A Docker stub run using `tests/fixtures/salmon_input.csv` completed through ISAR visualization, contrast summary, and MultiQC.
- The stub output correctly contained no `fastqc/`, `fastp/`, or `salmon/` directory, confirming that precomputed Salmon input bypasses read processing.

## Current Minimum Checks Before Future Commits

For normal implementation changes:

```bash
nextflow config -profile test,docker
nextflow run . -profile test,docker --outdir results/test
```

For changes touching GO enrichment workflow wiring:

```bash
nextflow run . \
    -profile test,docker \
    --run_go_enrichment \
    --outdir results/test_go_enrichment
```

For changes touching real GO enrichment behavior:

- Run a VM validation with a non-empty `isar/isar_analysis/go_gene_scores.csv`.
- Confirm `go/go_input/go_gene_scores.csv`, per-contrast gene lists, `go/go_enrichment.csv`, and `go/go_summary.txt` are produced.
- Confirm the GO enrichment container does not need `IsoformSwitchAnalyzeR`; it should consume the CSV handoff only.

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

For changes touching IUPred2A workflow wiring:

```bash
nextflow run . \
    -profile test,docker \
    -stub-run \
    --run_iupred2a \
    --outdir results/test_iupred2a_stub
```

For changes touching SignalP workflow wiring:

```bash
nextflow run . \
    -profile test,docker \
    -stub-run \
    --run_signalp \
    --signalp_container <licensed-signalp-5.0b-image-or-sif> \
    --outdir results/test_signalp_stub
```

For changes touching DeepTMHMM, DeepLoc2, or annotated switch plot workflow wiring:

```bash
nextflow run . \
    -profile test,docker \
    -stub-run \
    --run_deeptmhmm true \
    --run_deeploc2 true \
    --deeploc2_container <licensed-deeploc-2.1-image-or-sif> \
    --run_annotated_switch_plots true \
    --outdir results/test_milestone_b_stub
```

For changes touching the real SignalP command or import behavior:

- Run a VM validation on data that emits non-empty amino-acid FASTA input.
- The tiny `test,docker` dataset may skip the actual `SIGNALP_RUN` process when no candidate AA records are produced, so it is useful for graph health but not sufficient for command-level SignalP validation.

For changes touching SRA mode:

- Run at least one small real SRA validation on the VM.
- Confirm `prefetch`, `fasterq-dump`, compression, and downstream FASTQ processing complete.

## Known Testing Limitations

- Full real-data runs are manual and not part of CI.
- SRA tests depend on network access and public archive availability.
- Pfam tests require a prepared local Pfam database.
- Full SignalP validation depends on a user-provided, appropriately licensed SignalP 5.0b container.
- Full DeepTMHMM and DeepLoc2 validation is manual because the tools are too heavy for ordinary CI.
- The synthetic fixture is intentionally tiny and cannot prove biological correctness.
- More nf-test snapshots/assertions should be added once the output contract stabilizes.
