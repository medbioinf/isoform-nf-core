# Interface and Input Data Model Explained

Last updated: 2026-06-02

This document explains the user-facing interface of the pipeline: which files users provide, what each column means, how those inputs become Nextflow metadata, and how the parameters control the analysis.

## Why the Interface Matters

The interface is the contract between the user and the pipeline.

For a reusable pipeline, the interface should be:

- explicit enough that users know what to provide
- strict enough that mistakes are caught early
- flexible enough to work for different datasets
- stable enough that downstream code does not become fragile

This pipeline currently supports two input modes:

- FASTQ mode with `--input`
- SRA mode with `--sra_manifest`

Both modes are converted internally into the same data structure:

```nextflow
[ meta, reads ]
```

where `meta` is sample metadata and `reads` are one or two FASTQ files.

## Input Mode 1: FASTQ Samplesheet

FASTQ mode is for users who already have sequencing files on disk.

Command-line parameter:

```bash
--input samplesheet.csv
```

Example:

```csv
sample,condition,replicate,fastq_1,fastq_2,strandedness,batch
CONTROL_REP1,control,1,reads/control_R1.fastq.gz,reads/control_R2.fastq.gz,auto,batch1
TREATED_REP1,treated,1,reads/treated_R1.fastq.gz,reads/treated_R2.fastq.gz,auto,batch1
```

The schema is defined in:

```text
assets/schema_input.json
```

### `sample`

The `sample` column is the user-defined name of one biological sample.

Examples:

- `CONTROL_REP1`
- `PATIENT32_NORMAL`
- `PATIENT32_CANCER`

The same `sample` name can appear in multiple rows if the same biological sample was sequenced more than once, for example across several lanes.

Important: repeated rows with the same `sample` are treated as technical runs of the same biological sample, not as biological replicates.

### `condition`

The `condition` column describes the biological or experimental group.

Examples:

- `normal`
- `cancer`
- `control`
- `treated`
- `sepsis`

IsoformSwitchAnalyzeR compares conditions. In a simple two-condition design, it asks:

> Does isoform usage differ between condition A and condition B?

### `replicate`

The `replicate` column is the biological replicate number within a condition.

Examples:

- `1`
- `2`
- `3`

For statistical isoform switch testing, biological replicates are important because they let the statistical model estimate normal variation within each condition.

The current ISAR test path requires two conditions with at least two samples per condition. A `1 vs 1` comparison can be imported into ISAR, but the DEXSeq-based statistical test is skipped.

### `fastq_1`

The `fastq_1` column points to the first FASTQ file.

For single-end data, this is the only FASTQ file.

For paired-end data, this is read 1, often called `R1`.

### `fastq_2`

The `fastq_2` column points to the second FASTQ file for paired-end sequencing.

For single-end data, this column is empty.

The pipeline detects single-end vs paired-end from whether `fastq_2` is present.

### `strandedness`

The `strandedness` column tells the pipeline whether the RNA-seq library preserves information about which DNA strand the RNA came from.

Allowed values:

- `auto`
- `forward`
- `reverse`
- `unstranded`

Non-biologist explanation:

DNA has two strands. Genes can be encoded on either strand. Some RNA-seq protocols preserve the strand information, so we can tell which direction the original RNA came from. Other protocols lose this information.

Why it matters:

Salmon needs to know the library type to interpret reads correctly. Wrong strandedness can assign reads to the wrong transcript or reduce quantification quality.

Current behavior:

- `forward` maps to Salmon `SF` for single-end and `ISF` for paired-end.
- `reverse` maps to Salmon `SR` for single-end and `ISR` for paired-end.
- `unstranded` maps to Salmon `U` for single-end and `IU` for paired-end.
- `auto` currently behaves like the unstranded defaults unless `--salmon_lib_type` is set.

### `batch`

The `batch` column is optional metadata.

A batch is a technical grouping that might affect measurements, such as:

- sequencing date
- sequencing machine
- library preparation batch
- lab processing batch

The current V1 pipeline validates and carries the `batch` value, but the first ISAR implementation does not yet model batch effects statistically. It is included because it is useful metadata and gives us room to support batch-aware analysis later.

## Input Mode 2: SRA Manifest

SRA mode is for public sequencing data from SRA / ENA / DRA.

Command-line parameter:

```bash
--sra_manifest sra_manifest.csv
```

Example:

```csv
sample,condition,replicate,run_accession,strandedness,batch
PATIENT32_NORMAL,normal,1,SRR5275286,auto,gse95132
PATIENT32_CRC,crc,1,SRR5275287,auto,gse95132
```

The schema is defined in:

```text
assets/schema_sra_manifest.json
```

### Why SRA Mode Has `run_accession`

SRA stores public sequencing experiments as runs. A run accession such as `SRR5275286` identifies one sequencing run that can be downloaded.

The pipeline needs this accession because it directly calls:

```bash
prefetch SRR...
fasterq-dump SRR...
```

The SRA manifest uses the same biological metadata columns as FASTQ mode:

- `sample`
- `condition`
- `replicate`
- `strandedness`
- `batch`

The only difference is that instead of local FASTQ paths, SRA mode provides `run_accession`.

### Minimal vs Full SRA Mode

The current implementation is minimal SRA mode.

It expects users to provide run-level accessions such as:

- `SRR...`
- `ERR...`
- `DRR...`

It does not yet infer runs from higher-level accessions such as:

- `GSE...`
- `GSM...`
- `SRP...`

This is intentional for V1. Resolving higher-level accessions into runs requires additional metadata fetching and more edge-case handling.

## Contrast File

The contrast file is optional.

Command-line parameter:

```bash
--contrasts contrasts.csv
```

Example:

```csv
contrast,case,control
cancer_vs_normal,cancer,normal
```

The schema is defined in:

```text
assets/schema_contrasts.json
```

### `contrast`

The `contrast` column gives a human-readable name for the comparison.

Example:

```text
cancer_vs_normal
```

### `case`

The `case` column is the condition of interest.

Example:

```text
cancer
```

### `control`

The `control` column is the baseline condition.

Example:

```text
normal
```

Internally, the current R script converts this into the comparison format expected by IsoformSwitchAnalyzeR:

```r
condition_1 = control
condition_2 = case
```

## Reference Inputs

The pipeline needs reference files because V1 is reference-based.

### `--transcript_fasta`

This is a FASTA file containing transcript sequences.

Non-biologist explanation:

If a gene is like a recipe with optional sections, transcripts are the recipe variants. The transcript FASTA contains the actual sequence of each known recipe variant.

Used by:

- Salmon index creation
- IsoformSwitchAnalyzeR import

Important:

Transcript IDs should match between the transcript FASTA and the GTF file.

### `--gtf`

This is an annotation file describing where genes, transcripts, and exons are located.

Non-biologist explanation:

The GTF is like a map that says:

- this gene is on this chromosome
- this transcript belongs to this gene
- this transcript contains these exon pieces

Used by:

- Salmon `--geneMap`
- ISAR transcript-to-gene structure import
- filtering annotation to quantified transcripts

### `--genome_fasta`

This is optional.

It contains the genome sequence.

Used by:

- Salmon index creation as decoy sequence

Non-biologist explanation:

When quantifying transcripts, some reads may look like they match a transcript but actually come from other genomic regions. Providing the genome as decoy sequence helps Salmon avoid some incorrect transcript assignments.

The pipeline can run without `--genome_fasta`, but a decoy-aware index is usually more specific.

## Salmon Options

### `--salmon_lib_type`

This advanced parameter lets users directly specify Salmon's library type.

If unset, the pipeline derives the Salmon library type from:

- single-end vs paired-end
- `strandedness`

Examples of Salmon library types:

- `A`: automatic inference
- `U`: single-end unstranded
- `IU`: paired-end unstranded
- `SF`: single-end forward stranded
- `SR`: single-end reverse stranded
- `ISF`: paired-end inward-oriented forward stranded
- `ISR`: paired-end inward-oriented reverse stranded

For most users, the samplesheet `strandedness` column is easier than setting `--salmon_lib_type` manually.

## IsoformSwitchAnalyzeR Options

### `--run_isar`

Default:

```text
true
```

If true, the pipeline runs the ISAR import and statistical test after Salmon.

If false, the pipeline stops after Salmon and MultiQC.

This is useful for debugging quantification separately from downstream biological analysis.

### `--isar_dif_cutoff`

Default:

```text
0.1
```

`dIF` means difference in isoform fraction.

Isoform fraction is the fraction of a gene's expression that comes from one isoform.

Example:

- isoform uses 20 percent of gene expression in normal
- isoform uses 35 percent in cancer
- `dIF = 0.35 - 0.20 = 0.15`

With `--isar_dif_cutoff 0.1`, this isoform passes the effect-size threshold because `0.15` is larger than `0.1`.

### `--isar_qvalue_cutoff`

Default:

```text
0.05
```

The q-value is a multiple-testing corrected p-value.

Non-biologist explanation:

The pipeline tests many genes and isoforms. If we test thousands of things, some will look interesting by random chance. The q-value adjusts for this by estimating how many significant-looking results might be false positives.

A q-value cutoff of `0.05` means we accept candidates only if the corrected statistical support is strong enough under the chosen threshold.

### `--isar_top_n`

Default:

```text
10
```

Controls how many top switch candidates are written to `top_switches.csv`.

### `--run_isar_visualization`

Default:

```text
true
```

If true, the pipeline creates lightweight plots and candidate tables from the IsoformSwitchAnalyzeR result object. This step is designed to be safe for real smoke tests: if no statistically significant switches are detected, it writes explanatory notes and empty summary tables instead of crashing.

### `--isar_visualization_top_n`

Default:

```text
10
```

Controls how many top genes are used for visualization outputs such as top switching genes and per-gene isoform usage plots.

## MultiQC Options

### `--multiqc_title`

Sets the title shown in the MultiQC report.

This is useful when comparing many runs.

Example:

```bash
--multiqc_title "GSE95132 2v2 smoke test"
```

## How Input Rows Become Nextflow Metadata

Input parsing happens in:

```text
subworkflows/local/utils_nfcore_isoform_pipeline/main.nf
```

For FASTQ mode:

1. `samplesheetToList()` reads `params.input`.
2. `assets/schema_input.json` validates rows.
3. `sampleRowToInput()` converts each row into an internal tuple.
4. `groupTuple()` groups repeated rows with the same sample name.
5. `validateInputSamplesheet()` checks repeated runs are consistent.
6. The workflow emits `ch_samplesheet`.

For SRA mode:

1. `samplesheetToList()` reads `params.sra_manifest`.
2. `assets/schema_sra_manifest.json` validates rows.
3. `sraRowToInput()` converts each row into a run tuple.
4. `groupTuple()` groups repeated runs by sample.
5. `validateInputSamplesheet()` checks repeated runs are consistent.
6. `flatMap()` expands runs again so each SRA accession can be downloaded separately.
7. The workflow emits `ch_sra_manifest`.

## Validation Rules

The pipeline fails early if:

- both `--input` and `--sra_manifest` are provided
- neither `--input` nor `--sra_manifest` is provided
- `--transcript_fasta` is missing
- `--gtf` is missing
- repeated runs of the same sample mix single-end and paired-end data
- repeated runs of the same sample have conflicting metadata
- a contrast references a condition not present in the input

This early validation is important because otherwise errors would appear much later inside Salmon, ISAR, or a container task.

## Current Interface Limitations

The current V1 interface intentionally does not support every possible RNA-seq design.

Known limitations:

- The statistical ISAR test path is currently designed for two-condition comparisons.
- Batch metadata is carried and validated but not yet modeled.
- SRA mode requires run accessions and does not resolve `GSE` or `GSM` accessions automatically.
- Novel isoform discovery is not part of V1.
- Long-read RNA-seq input is not part of V1.
- Visualization options from the prototype pipeline are not yet exposed.

These limitations are not necessarily design flaws. They keep the first reusable version understandable and testable.
