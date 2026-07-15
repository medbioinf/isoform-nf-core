# Tools and Biology Explained

Last updated: 2026-07-05

This document explains the biological ideas and tools used by the current pipeline. It is written for readers who are not biologists.

## The Biological Problem

Human genes are not always used in only one fixed way.

A gene can produce multiple transcript variants, called isoforms. These isoforms can include or skip different exons, use different start or end sites, or produce proteins with different functional regions.

An isoform switch happens when two conditions use different isoforms of the same gene.

Example:

- In normal tissue, gene X mostly uses isoform A.
- In cancer tissue, gene X mostly uses isoform B.
- The total amount of gene X might stay similar, but the product changes.

That can matter because different isoforms may produce proteins with different domains, localization, membrane topology, disorder, or stability.

## Gene vs Transcript vs Isoform

### Gene

A gene is a region of DNA that can produce RNA.

Non-biologist analogy:

Think of a gene as a recipe with optional paragraphs.

### Transcript

A transcript is an RNA molecule produced from a gene.

Non-biologist analogy:

A transcript is one concrete version of the recipe after choosing which optional paragraphs to include.

### Isoform

In this project, isoform and transcript are often used almost interchangeably. More precisely, an isoform is a transcript variant produced from the same gene.

Non-biologist analogy:

If the gene is the base recipe, isoforms are recipe variants.

## Exons, Introns, and Splicing

Genes are often split into exons and introns.

- Exons are sequence segments that can be retained in mature RNA.
- Introns are removed during RNA processing.
- Splicing is the process that removes introns and joins exons.

Alternative splicing means different exon combinations can be produced from the same gene.

This is one major reason genes can have multiple isoforms.

## Short Reads and Paired-end Reads

The current pipeline focuses on short-read RNA-seq, usually Illumina-style data.

### Single-end

Single-end sequencing reads one end of each RNA-derived fragment.

Pipeline representation:

```csv
fastq_1,fastq_2
sample_R1.fastq.gz,
```

### Paired-end

Paired-end sequencing reads both ends of each fragment.

Pipeline representation:

```csv
fastq_1,fastq_2
sample_R1.fastq.gz,sample_R2.fastq.gz
```

Paired-end data usually gives more information about where a fragment came from.

## Reference-based Analysis

The current pipeline is reference-based.

That means it uses known transcript sequences and annotations:

- transcript FASTA
- GTF annotation
- optional genome FASTA

It does not currently discover completely new isoforms that are absent from the reference.

Why:

Salmon quantifies abundance against a known transcriptome. If a transcript is not in the reference transcript FASTA, Salmon cannot directly quantify it as its own transcript.

Novel isoform detection would require a different path, such as transcript reconstruction or long-read transcript discovery, before quantification and ISAR analysis.

## Tool 1: SRA Tools

Used in:

```text
modules/local/fetch_sra_fastq
```

Tools:

- `prefetch`
- `fasterq-dump`

Purpose:

Download public sequencing runs and convert them to FASTQ.

Non-biologist explanation:

SRA is a public archive for sequencing data. It stores data in archive formats, not always directly as the FASTQ files most tools expect. SRA tools convert public run accessions into usable FASTQ files.

Current pipeline behavior:

1. `prefetch` downloads the run archive.
2. `fasterq-dump` converts the archive to FASTQ.
3. `pigz` compresses the FASTQ files.
4. The pipeline continues as if the user had supplied FASTQ files directly.

## Tool 2: FastQC

Used in:

```text
modules/nf-core/fastqc
```

Purpose:

Quality control of raw sequencing reads.

Non-biologist explanation:

FastQC is like an inspection report for raw sequencing files. It does not fix the data. It tells us whether the raw data looks healthy or suspicious.

FastQC checks things such as:

- base quality scores across the read
- GC content
- adapter contamination
- overrepresented sequences
- read length distribution
- duplicated sequences

What users do with the information:

- Check whether sequencing data is usable.
- Notice strong adapter contamination.
- Notice poor quality at read ends.
- Compare samples for suspicious outliers.
- Include QC evidence in reports.

In this pipeline, FastQC reports are collected by MultiQC.

## Tool 3: cat/fastq

Used in:

```text
modules/nf-core/cat/fastq
```

Purpose:

Concatenate multiple FASTQ runs belonging to the same biological sample.

Non-biologist explanation:

Sometimes the same sample is sequenced in multiple lanes or runs. That creates multiple FASTQ files for one biological sample. For downstream analysis, these should usually be combined so the sample is counted once.

Important distinction:

- Multiple runs of the same `sample` are technical repeats and are concatenated.
- Different biological replicates must have different `sample` names.

## Tool 4: fastp

Used in:

```text
modules/nf-core/fastp
```

Purpose:

Trim adapters and low-quality sequence before quantification.

Non-biologist explanation:

Sequencing reads can contain unwanted technical sequence or low-confidence bases. fastp cleans the reads so quantification is less affected by technical noise.

fastp can:

- detect adapter contamination
- trim adapters
- remove or trim low-quality bases
- produce HTML and JSON QC reports
- detect adapters for paired-end reads

The cleaned reads are passed to Salmon.

## Tool 5: Salmon

Used in:

```text
modules/nf-core/salmon/index
modules/nf-core/salmon/quant
```

Purpose:

Estimate transcript abundance from RNA-seq reads.

Non-biologist explanation:

Salmon asks:

> Given these sequencing reads and this catalog of known transcript sequences, how much of each transcript is likely present in each sample?

It produces `quant.sf`, which contains transcript-level abundance estimates.

Important columns in `quant.sf`:

- `Name`: transcript ID.
- `Length`: transcript length.
- `EffectiveLength`: adjusted transcript length used for quantification.
- `TPM`: normalized transcript abundance.
- `NumReads`: estimated number of reads assigned to the transcript.

### Why Salmon Needs a Transcript FASTA

The transcript FASTA is the catalog of known transcript sequences.

Salmon builds an index from it, similar to making a searchable database.

### Why Salmon Uses the GTF

The GTF maps transcripts back to genes.

This allows Salmon to also produce gene-level summaries and helps downstream tools understand which isoforms belong to which gene.

### Why the Genome FASTA Is Optional

If provided, the genome FASTA is used as decoy sequence during indexing.

This helps Salmon avoid incorrectly assigning reads to transcripts when they better match non-transcript genomic sequence.

## Tool 6: IsoformSwitchAnalyzeR

Used in:

```text
modules/local/isar_analysis
bin/run_isar_analysis.R
```

Purpose:

Analyze whether isoform usage changes between conditions.

Non-biologist explanation:

Salmon tells us how much each transcript appears in each sample. IsoformSwitchAnalyzeR looks at genes with multiple isoforms and asks whether the relative usage of those isoforms changes between conditions.

Example:

```text
Gene Z
Normal: 80 percent isoform A, 20 percent isoform B
Cancer: 45 percent isoform A, 55 percent isoform B
```

This may be an isoform switch.

### Isoform Fraction

Isoform fraction is the share of a gene's expression that comes from one isoform.

Example:

If a gene has two isoforms:

- isoform A expression = 80
- isoform B expression = 20
- total gene expression = 100

Then:

- isoform A fraction = 0.8
- isoform B fraction = 0.2

### dIF

`dIF` means difference in isoform fraction.

Example:

- isoform fraction in normal = 0.20
- isoform fraction in cancer = 0.45
- `dIF = 0.45 - 0.20 = 0.25`

This means the isoform is used 25 percentage points more in cancer than normal.

### q-value

The q-value is a multiple-testing corrected p-value.

Non-biologist explanation:

If we test thousands of isoforms, some will look interesting by chance. The q-value corrects for the fact that many tests were performed.

Lower q-values mean stronger statistical support.

### Cutoffs and Top Candidates

The pipeline first applies both switch cutoffs:

```text
q-value < q-value cutoff
absolute dIF >= dIF cutoff
```

Only isoforms that pass both filters are treated as significant switch candidates.

Top-N outputs are then ranked from this filtered candidate set. In most candidate tables, the strongest candidates are sorted by:

1. lowest isoform switch q-value
2. largest absolute `dIF`

This means an isoform with a very strong q-value but a tiny `dIF` is excluded before ranking if it does not pass the `dIF` cutoff.

### Why Replicates Matter

IsoformSwitchAnalyzeR needs biological replicates to estimate normal variation inside each condition.

With only one sample per condition, it is impossible to tell whether an observed difference is a real condition effect or just random sample-to-sample variation.

That is why the current pipeline imports `1 vs 1` data but skips DEXSeq-based statistical testing.

### DEXSeq

DEXSeq is the statistical method currently used by IsoformSwitchAnalyzeR for testing differential isoform usage.

Non-biologist explanation:

DEXSeq helps answer:

> Is the isoform usage difference between conditions larger than expected from normal replicate variation?

## Tool 7: WebGestaltR

Used in:

```text
modules/local/isoform_go_enrichment
bin/run_go_enrichment.R
```

Purpose:

Run gene ontology or other gene-set enrichment analyses for genes with significant isoform switches.

Non-biologist explanation:

IsoformSwitchAnalyzeR finds individual switching genes and isoforms. WebGestaltR asks whether those switching genes share broader biological themes, such as immune response, cell adhesion, or membrane organization.

Important distinction:

GO enrichment is gene-level. It does not explain which exact isoform gained or lost a protein domain. It complements the annotated switch plots, which remain the better view for detailed isoform-level consequences.

## Tool 8: MultiQC

Used in:

```text
modules/nf-core/multiqc
```

Purpose:

Collect many tool reports into one HTML report.

Non-biologist explanation:

Instead of opening one FastQC file per sample, one fastp report per sample, and Salmon logs separately, MultiQC summarizes the run in one report.

It is mainly a reporting and QC aggregation tool.

## Optional Annotation Tools

These tools are integrated as optional modules. They run after IsoformSwitchAnalyzeR and add biological interpretation to significant isoform switch candidates. Their outputs are available as standalone tables and summary plots, and the imported annotations can also be used by the annotated switch plot module.

### Pfam

Purpose:

Detect protein domains.

Non-biologist explanation:

Protein domains are functional building blocks. If an isoform switch removes or adds a domain, the resulting protein may behave differently.

Why useful:

Pfam can help explain the biological consequence of an isoform switch.

Operational note:

Pfam needs a large database and additional runtime handling.

### SignalP

Purpose:

Predict signal peptides.

Non-biologist explanation:

A signal peptide is like an address label that sends a protein into the secretory pathway. The prediction does not by itself determine the protein's final location or prove that it is secreted.

Why useful:

An isoform switch could add or remove this address label.

Operational note:

SignalP 5.0b must be obtained under the appropriate DTU license and supplied in a private/local container; the pipeline does not redistribute it.

### IUPred2A

Purpose:

Predict intrinsically disordered protein regions.

Non-biologist explanation:

Some protein regions do not fold into one stable shape. These flexible regions can matter for regulation and interactions.

Why useful:

An isoform switch could add or remove disordered regions.

Operational note:

The current container is usable for the pipeline, but long-term public release should still review container provenance and runtime reporting behavior.

### DeepTMHMM

Purpose:

Predict transmembrane topology.

Non-biologist explanation:

Some proteins sit in membranes. DeepTMHMM predicts which parts cross the membrane and which parts are inside or outside the cell or organelle.

Why useful:

An isoform switch could change whether a protein is membrane-bound or which side of a membrane a region faces.

Operational note:

DeepTMHMM can be slow on real datasets because model inference is more expensive than simple table conversion.

### DeepLoc2

Purpose:

Predict subcellular localization.

Non-biologist explanation:

Predicts where a protein likely goes inside the cell, such as nucleus, cytoplasm, mitochondria, or membrane.

Why useful:

An isoform switch could change where the protein acts.

Operational note:

DeepLoc2 may need model assets and can require more runtime on first use depending on the container cache state.

## What Results Mean in Practice

A significant isoform switch candidate usually means:

- a gene has at least two isoforms
- the relative usage of those isoforms differs between conditions
- the difference is large enough to pass the `dIF` cutoff
- the statistical evidence is strong enough to pass the q-value cutoff

It does not automatically mean:

- the switch is biologically causal
- the protein definitely changes function
- the result is valid without checking QC
- the finding is clinically meaningful

The pipeline produces candidates. Biological interpretation still needs domain knowledge, literature review, visualization, and ideally validation.
