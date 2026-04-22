# V1 Interface Design

This document defines the proposed user-facing interface for the first real version of the `isoform-nf-core` pipeline.

It is intended to be the contract that guides the first implementation steps:

- `nextflow.config` parameter design
- `nextflow_schema.json`
- input validation
- workflow branching
- test data and examples

## V1 scope

The first version of the pipeline should support the following analysis path:

`FASTQ or optional SRA manifest -> QC / trimming -> Salmon -> IsoformSwitchAnalyzeR`

V1 should support:

- paired-end RNA-seq
- single-end RNA-seq
- primary input via local FASTQ samplesheet
- optional input via SRA manifest
- transcript-level quantification with Salmon
- cohort-level isoform-switch analysis with `IsoformSwitchAnalyzeR`

V1 should not include:

- dataset-branded modes such as `gse95132_crc`
- bundled GEO-specific manifests as defaults
- hardcoded condition labels such as `normal` and `crc`
- every optional external ISAR annotation tool from the vignette

## Input modes

The pipeline should support exactly two top-level input modes.

### Mode 1: FASTQ samplesheet

This should be the primary and recommended mode.

User provides:

- `--input`
- `--transcript_fasta`
- `--genome_fasta`
- `--gtf`

The samplesheet should describe the cohort and FASTQ locations directly.

### Mode 2: SRA manifest

This should be an optional convenience mode for public datasets.

User provides:

- `--sra_manifest`
- `--transcript_fasta`
- `--genome_fasta`
- `--gtf`

The manifest should describe the cohort and the SRA accessions.

The pipeline should then:

1. download / convert public data to FASTQ
2. normalize those inputs into the same internal representation used by FASTQ mode
3. continue through the same generic workflow

## Input mode rules

The pipeline should enforce:

- exactly one of `--input` or `--sra_manifest`
- reference files are always required in v1
- the rest of the workflow should operate on one normalized internal sample representation

## FASTQ samplesheet design

### Required columns

- `sample`
- `condition`
- `fastq_1`

### Optional columns

- `fastq_2`
- `layout`
- `patient_id`

### Column meaning

- `sample`
  - unique sample identifier used throughout the pipeline
- `condition`
  - biological condition used for comparison in the ISAR step
- `fastq_1`
  - path to the first FASTQ file, or the only FASTQ file for single-end data
- `fastq_2`
  - path to the second FASTQ file for paired-end data
- `layout`
  - optional explicit layout, one of `single` or `paired`
- `patient_id`
  - optional blocking / pairing variable for matched designs

### Validation rules

- `sample` must be unique
- `condition` must be non-empty
- `fastq_1` must always be present
- if `layout == paired`, then `fastq_2` must be present
- if `layout` is missing:
  - infer `paired` when `fastq_2` is populated
  - otherwise infer `single`

### Example FASTQ samplesheet

```csv
sample,condition,fastq_1,fastq_2,layout,patient_id
sample_normal_1,normal,/data/sample_normal_1_R1.fastq.gz,/data/sample_normal_1_R2.fastq.gz,paired,patient01
sample_tumor_1,tumor,/data/sample_tumor_1_R1.fastq.gz,/data/sample_tumor_1_R2.fastq.gz,paired,patient01
sample_crc_2,crc,/data/sample_crc_2.fastq.gz,,single,patient02
```

## SRA manifest design

### Required columns

- `sample`
- `condition`
- `sra_run`

### Optional columns

- `layout`
- `patient_id`

### Column meaning

- `sample`
  - unique sample identifier used throughout the pipeline
- `condition`
  - biological condition used for comparison in the ISAR step
- `sra_run`
  - SRA run accession such as `SRR1234567`
- `layout`
  - optional explicit layout, one of `single` or `paired`
- `patient_id`
  - optional blocking / pairing variable

### Validation rules

- `sample` must be unique
- `condition` must be non-empty
- `sra_run` must be non-empty
- if `layout` is supplied, it must be `single` or `paired`
- if `layout` is missing, allow the downloader / conversion layer to determine the actual output shape

### Example SRA manifest

```csv
sample,condition,sra_run,layout,patient_id
patient32_normal,normal,SRR5275286,single,32
patient32_crc,crc,SRR5275287,single,32
```

## Reference inputs

V1 should require these references in both modes:

- `--transcript_fasta`
- `--genome_fasta`
- `--gtf`

### Why all three are required

- `transcript_fasta` and `genome_fasta` are needed for decoy-aware Salmon indexing
- `gtf` is needed for the ISAR import / annotation layer

The pipeline should not hide dataset-specific reference choices behind branded flags.

## Proposed initial user-facing parameters

These are the parameters that should likely appear in the first real schema revision.

### Input and output

- `input`
  - FASTQ samplesheet path
- `sra_manifest`
  - SRA manifest path
- `outdir`
  - output directory

### References

- `transcript_fasta`
- `genome_fasta`
- `gtf`

### Workflow behavior

- `skip_qc`
  - optional boolean to skip `FASTQC`
- `skip_trimming`
  - optional boolean to skip `FASTP`
- `run_isar`
  - boolean, default `true`

### Salmon options

- `salmon_index_k`
  - k-mer size used when building the Salmon index
  - smaller values can help with shorter reads
- `salmon_quant_libtype`
  - default likely `A`
  - tells Salmon what library type to assume, or to auto-detect it
- `salmon_fld_mean`
  - relevant for single-end mode
  - expected average fragment length for single-end data
- `salmon_fld_sd`
  - relevant for single-end mode
  - expected fragment length spread for single-end data
- `salmon_seq_bias`
  - whether to enable sequence-specific bias correction
- `salmon_gc_bias`
  - whether to enable GC-bias correction
- `salmon_pos_bias`
  - whether to enable positional bias correction

### ISAR comparison options

- `isar_condition_1`
  - first condition label to compare, for example `normal`
- `isar_condition_2`
  - second condition label to compare, for example `tumor`
- `dif_cutoff`
  - minimum change in isoform usage required before a switch is considered meaningful
- `qvalue_cutoff`
  - false-discovery-rate threshold used to decide statistical significance
- `isar_top_n`
  - how many top switch candidates to include in the main ranked output
- `isar_detect_unwanted_effects`
  - whether ISAR should try to detect unwanted technical effects during import

### SRA mode options

- `sra_prefetch_max_size`
  - maximum archive size that `prefetch` is allowed to download
- `download_threads`
  - number of threads to use for download and FASTQ conversion steps

### Reporting

- `multiqc_title`
  - custom title shown at the top of the MultiQC report
  
## Internal normalized sample model

Regardless of input mode, the pipeline should convert inputs into one shared internal structure.

Suggested meta fields:

- `id`
- `condition`
- `layout`
- `patient_id`
- `sra_run`

Suggested read representation:

- single-end: one FASTQ file
- paired-end: ordered pair of FASTQ files

This normalization is important because it lets the main workflow stay generic.

## Expected V1 outputs

The pipeline should produce at least:

- per-sample `FASTQC` results unless skipped
- per-sample `FASTP` reports unless skipped
- `MultiQC` report
- per-sample Salmon quantification directories
- cohort metadata / design outputs used by ISAR
- ISAR result tables
- ISAR plots
- pipeline execution metadata

## Immediate implementation consequences

Before porting large amounts of code, the next steps should be:

1. update `nextflow.config` params in `isoform-nf-core`
2. replace the template `input` schema with the proposed FASTQ samplesheet contract
3. add a second schema entry for `sra_manifest`
4. define the internal sample metadata contract in the workflow code
5. only then start porting reusable logic from the prototype repo
