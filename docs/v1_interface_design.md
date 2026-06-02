# V1 Interface Status

Last updated: 2026-06-02

This document records the current implemented V1 interface. It replaces the earlier planning notes and should be read together with `docs/usage.md` and `nextflow_schema.json`.

## Scope

The implemented V1 path is:

`FASTQ or SRA manifest -> FastQC -> cat/fastq -> fastp -> Salmon -> IsoformSwitchAnalyzeR -> MultiQC`

The pipeline is reference-based. It quantifies transcripts from the supplied reference transcript FASTA and GTF; it does not discover de novo novel isoforms.

## Input Modes

Exactly one input mode must be provided:

- `--input`: local FASTQ samplesheet
- `--sra_manifest`: public SRA/ENA/DRA run manifest

Both modes are normalized into the same downstream read-processing path.

## FASTQ Samplesheet

Required columns:

- `sample`
- `condition`
- `replicate`
- `fastq_1`
- `strandedness`

Optional columns:

- `fastq_2`
- `batch`

Example:

```csv
sample,condition,replicate,fastq_1,fastq_2,strandedness,batch
CONTROL_REP1,control,1,/data/control_R1.fastq.gz,/data/control_R2.fastq.gz,auto,batch1
TREATED_REP1,treated,1,/data/treated_R1.fastq.gz,,auto,batch1
```

Rows with the same `sample` are treated as multiple sequencing runs for one biological sample. They must agree on `condition`, `replicate`, `strandedness`, `batch`, and read layout.

## SRA Manifest

Required columns:

- `sample`
- `condition`
- `replicate`
- `run_accession`
- `strandedness`

Optional columns:

- `batch`

Example:

```csv
sample,condition,replicate,run_accession,strandedness,batch
CONTROL_REP1,control,1,SRR000001,auto,batch1
TREATED_REP1,treated,1,SRR000002,auto,batch1
```

SRA mode uses `prefetch` followed by `fasterq-dump`, compresses the generated FASTQs, detects single-end versus paired-end output, and then enters the normal FASTQ path.

## Contrast File

The optional `--contrasts` file defines pairwise condition comparisons:

```csv
contrast,case,control
treated_vs_control,treated,control
```

The `case` and `control` values must exist in the input `condition` column.

## Reference Inputs

Required:

- `--transcript_fasta`
- `--gtf`

Optional:

- `--genome_fasta`

If `--genome_fasta` is provided, Salmon builds a decoy-aware index. Transcript IDs should match between the transcript FASTA and GTF.

## Main Analysis Parameters

- `--salmon_lib_type`: optional explicit Salmon library type; leave unset for pipeline-derived/default behavior.
- `--run_isar`: run or skip IsoformSwitchAnalyzeR, default `true`.
- `--isar_dif_cutoff`: minimum absolute isoform fraction difference for switch testing, default `0.1`.
- `--isar_qvalue_cutoff`: q-value cutoff for switch testing, default `0.05`.
- `--isar_top_n`: number of top switches to export when switches are detected.
- `--sra_prefetch_max_size`: maximum archive size accepted by SRA `prefetch`, default `100G`.

## Current Limitations

- `batch` is validated and preserved as metadata, but not yet modeled statistically in ISAR.
- The current ISAR wrapper runs the DEXSeq-based switch test only when there are exactly two conditions with at least two samples per condition.
- Additional IsoformSwitchAnalyzeR companion tools such as Pfam, SignalP, IUPred2A, DeepLoc2, and DeepTMHMM are not part of the reusable V1 path yet.
