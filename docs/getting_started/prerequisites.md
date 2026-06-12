# Prerequisites to Run the Pipeline

Last updated: 2026-06-12

This document explains what a user needs before running the pipeline. It focuses on practical setup: software runtime, input data, reference files, optional Pfam annotation resources, and compute/storage expectations.

## Short Version

For the normal reproducible path, the user should bring:

- Nextflow plus Java.
- A container runtime such as Docker, Apptainer, or Singularity.
- Either a FASTQ samplesheet or an SRA manifest.
- A transcript FASTA.
- A GTF annotation file.
- Optionally, a genome FASTA for Salmon decoy-aware indexing.
- Optionally, a contrast file if the dataset has more than one comparison.
- Optionally, a prepared Pfam database directory if Pfam domain annotation should be run.

The user does **not** need to manually install FastQC, fastp, Salmon, SRA Toolkit, MultiQC, Pfam scan, or IsoformSwitchAnalyzeR when using a container profile. Those tools run inside containers selected by the pipeline.

## Runtime Prerequisites

### Nextflow

The pipeline is a Nextflow workflow, so the machine that launches it needs Nextflow installed.

Check with:

```bash
nextflow -version
```

Nextflow itself needs Java. If `nextflow -version` fails with a Java error, install a Java runtime first.

Check Java with:

```bash
java -version
```

### Container Runtime

Use a container profile for reproducibility.

Recommended local or VM profile:

```bash
-profile docker
```

Recommended HPC profiles usually use one of:

```bash
-profile apptainer
-profile singularity
```

Why this matters:

- Containers define the tool versions used by each step.
- Users do not have to install bioinformatics tools manually.
- Runs are easier to reproduce on another machine.

If no container profile is used, Nextflow expects all tools to already exist on the system `PATH`. That is not recommended for this project.

### Internet Access

Internet access may be needed for:

- downloading the pipeline code from GitHub
- downloading containers on the first run
- downloading SRA data when using `--sra_manifest`

After containers and public data are cached, some runs can be repeated with less network access. On HPC systems, it is common to pre-download containers and references before running the workflow.

### Git

Git is useful when running directly from GitHub or when developing the pipeline locally.

Users who run a released pipeline through Nextflow usually do not interact with Git directly, but developers will.

## Input Data Prerequisites

The pipeline currently supports two input modes.

### FASTQ Mode

Use FASTQ mode when sequencing files already exist locally.

Required parameter:

```bash
--input samplesheet.csv
```

The samplesheet points to `.fastq.gz` or `.fq.gz` files and provides biological metadata such as sample name, condition, replicate, and strandedness.

Important expectations:

- FASTQ files should be gzip-compressed.
- Paired-end data uses both `fastq_1` and `fastq_2`.
- Single-end data leaves `fastq_2` empty.
- Multiple rows with the same `sample` are treated as technical runs of the same sample and are concatenated.

### SRA Mode

Use SRA mode when starting from public run accessions.

Required parameter:

```bash
--sra_manifest sra_manifest.csv
```

The SRA manifest provides run accessions such as `SRR...`, `ERR...`, or `DRR...`.

The pipeline downloads each run with SRA Toolkit inside a container, converts it to FASTQ, compresses the FASTQ files, and then continues through the same workflow path as FASTQ mode.

Current limitation:

- The pipeline expects run accessions.
- It does not yet resolve higher-level accessions such as `GSE...`, `GSM...`, or `SRP...` automatically.

## Reference File Prerequisites

The reference files define which transcripts can be quantified and analyzed.

### Transcript FASTA

Required parameter:

```bash
--transcript_fasta transcripts.fa.gz
```

This file contains the transcript sequences. Salmon uses it to build the transcriptome index and quantify transcript abundance.

Important:

- The transcript IDs should match the IDs in the GTF.
- Reference-based quantification can only quantify transcripts present in this reference.
- This pipeline version does not discover novel isoforms de novo.

### GTF Annotation

Required parameter:

```bash
--gtf annotation.gtf.gz
```

The GTF describes genes, transcripts, exons, and their relationships.

The pipeline uses it to connect transcript-level Salmon results back to genes and to prepare IsoformSwitchAnalyzeR input.

Important:

- The transcript IDs in the GTF should match the transcript FASTA.
- The gene/transcript annotation should come from the same reference release as the transcript FASTA.

### Genome FASTA

Optional but recommended parameter:

```bash
--genome_fasta genome.fa.gz
```

The genome FASTA can be used as decoy sequence when building the Salmon index. In simple terms, it helps Salmon distinguish true transcript matches from reads that map better to other genomic sequence.

For minimal tests, the pipeline can run without it. For real analyses, providing a matching genome FASTA is recommended.

## Experimental Design Prerequisites

### Conditions and Replicates

IsoformSwitchAnalyzeR needs biological conditions to compare.

For statistical switch testing, the current path expects:

- two or more conditions
- at least two biological samples per condition for the tested comparison

Very small smoke-test datasets can still be imported, but statistical testing may be skipped.

### Contrast File

Optional parameter:

```bash
--contrasts contrasts.csv
```

A contrast file is useful when the dataset contains more than one condition or when the desired comparison should be explicit.

Example:

```csv
contrast,case,control
cancer_vs_normal,cancer,normal
treated_vs_control,treated,control
```

If there are exactly two conditions, a contrast file may be optional. For clearer and more reproducible analyses, it is still a good idea to provide one.

## Optional Pfam Prerequisites

Pfam annotation is optional. It adds protein-domain information to significant isoform switches.

Normal user-facing parameter:

```bash
--pfam_db /path/to/pfam/Pfam37.0
```

The directory should contain a prepared Pfam database, including:

```text
Pfam-A.hmm
Pfam-A.hmm.dat
Pfam-A.hmm.h3f
Pfam-A.hmm.h3i
Pfam-A.hmm.h3m
Pfam-A.hmm.h3p
```

The `.h3*` files are HMMER index files generated with:

```bash
hmmpress Pfam-A.hmm
```

If index files are missing, the pipeline currently tries to run `hmmpress`. For shared HPC/VM use, it is better to prepare the Pfam database once and reuse it.

Why the database is not downloaded by default:

- Pfam files are large.
- Downloads can be slow or unavailable on clusters.
- Automatically using the latest database can hurt reproducibility.
- A prepared database can be reused across many analyses.

There is also an advanced `--pfam_results` option for importing an existing `pfam_scan.pl` result file. This is mainly useful for debugging, development, or rerunning downstream plots without repeating the expensive Pfam scan. Normal users should prefer `--pfam_db`.

## Compute and Storage Expectations

The pipeline can run small tests on a laptop, but real RNA-seq datasets are better run on a VM or HPC system.

Expect disk usage from:

- downloaded SRA archives
- converted FASTQ files
- Nextflow `work/` directories
- Salmon indices and quantification outputs
- ISAR result objects
- optional Pfam scan outputs

Pfam scanning can be noticeably slower than lightweight plotting steps. In our GSE50760 validation subset, the real Pfam scan and import were substantial enough to treat Pfam as an annotation module, not just a tiny visualization add-on.

Practical advice:

- Run real datasets on a VM/HPC rather than a laptop.
- Keep `work/` on a disk with enough space.
- Use `-resume` when restarting interrupted runs.
- Keep reference files and Pfam databases in stable shared locations.

## Minimal Example Commands

FASTQ input:

```bash
nextflow run Anton-Bch/isoform-nf-core \
    --input samplesheet.csv \
    --contrasts contrasts.csv \
    --transcript_fasta reference/transcripts.fa.gz \
    --genome_fasta reference/genome.fa.gz \
    --gtf reference/annotation.gtf.gz \
    --outdir results \
    -profile docker
```

SRA input:

```bash
nextflow run Anton-Bch/isoform-nf-core \
    --sra_manifest sra_manifest.csv \
    --contrasts contrasts.csv \
    --transcript_fasta reference/transcripts.fa.gz \
    --genome_fasta reference/genome.fa.gz \
    --gtf reference/annotation.gtf.gz \
    --outdir results \
    -profile docker
```

With optional Pfam annotation:

```bash
nextflow run Anton-Bch/isoform-nf-core \
    --input samplesheet.csv \
    --contrasts contrasts.csv \
    --transcript_fasta reference/transcripts.fa.gz \
    --genome_fasta reference/genome.fa.gz \
    --gtf reference/annotation.gtf.gz \
    --pfam_db reference/pfam/Pfam37.0 \
    --outdir results \
    -profile docker
```

## What Users Do Not Need to Install Manually

When using a container profile, users do not need local installations of:

- FastQC
- fastp
- Salmon
- SRA Toolkit
- MultiQC
- IsoformSwitchAnalyzeR
- Pfam scan
- HMMER

The user provides data and references. The pipeline provides the workflow logic and launches the required tools in containers.
