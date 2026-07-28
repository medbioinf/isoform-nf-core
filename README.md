# Anton-Bch/isoform-nf-core

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/Anton-Bch/isoform-nf-core)
[![GitHub Actions CI Status](https://github.com/Anton-Bch/isoform-nf-core/actions/workflows/nf-test.yml/badge.svg)](https://github.com/Anton-Bch/isoform-nf-core/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/Anton-Bch/isoform-nf-core/actions/workflows/linting.yml/badge.svg)](https://github.com/Anton-Bch/isoform-nf-core/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/Anton-Bch/isoform-nf-core)

## Introduction

**Anton-Bch/isoform-nf-core** is a reference-based RNA-seq pipeline for isoform switch analysis. It accepts local FASTQ samplesheets, public SRA run manifests, or existing Salmon quantification directories. Starting from reads, it performs QC and preprocessing, quantifies transcript abundance with Salmon, and imports the results into IsoformSwitchAnalyzeR for differential isoform usage analysis. Existing Salmon results enter directly at the IsoformSwitchAnalyzeR stage.

The current workflow runs:

1. Raw read QC with [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Optional SRA download and conversion with `prefetch` / `fasterq-dump`
3. Multiple-run concatenation with `cat/fastq`
4. Adapter and quality filtering with [`fastp`](https://github.com/OpenGene/fastp)
5. Transcript quantification with [`Salmon`](https://salmon.readthedocs.io/)
6. Isoform switch import and testing with [`IsoformSwitchAnalyzeR`](https://bioconductor.org/packages/IsoformSwitchAnalyzeR/)
7. ISAR summary visualizations, including per-comparison counts and UpSet-style intersections
8. Optional GO enrichment for significant switching genes with WebGestaltR
9. Optional functional annotation with Pfam, IUPred2A, SignalP, DeepTMHMM, and DeepLoc2
10. Optional annotated gene-level switch plots
11. Aggregated reporting with [`MultiQC`](http://multiqc.info/)

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,condition,replicate,fastq_1,fastq_2,strandedness,batch
CONTROL_REP1,control,1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,auto,batch1
TREATED_REP1,treated,1,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,auto,batch1
```

Each row represents a FastQ file (single-end) or a pair of FastQ files (paired-end), plus the sample metadata needed to run different datasets without changing the pipeline code.

Alternatively, provide `--sra_manifest` with SRA run accessions and the same sample metadata columns. Minimal SRA mode downloads each run with `prefetch`, converts it with `fasterq-dump`, and then continues through the normal FASTQ path.

If compatible Salmon `quant.sf` results already exist, provide `--salmon_input` with sample metadata and one quantification directory per sample. This mode skips FastQC, fastp, Salmon indexing, and Salmon quantification. See [Interface and input data model explained](docs/getting_started/interface_and_inputs_explained.md) for the three input formats and their requirements.

Optionally, prepare a contrast file with pairwise condition comparisons:

`contrasts.csv`:

```csv
contrast,case,control
treated_vs_control,treated,control
```

Now, you can run the pipeline using:

```bash
nextflow run Anton-Bch/isoform-nf-core \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --contrasts contrasts.csv \
   --transcript_fasta reference/transcripts.fa.gz \
   --genome_fasta reference/genome.fa.gz \
   --gtf reference/annotation.gtf.gz \
   --outdir <OUTDIR>
```

The core workflow runs read QC, optional run concatenation, fastp preprocessing, Salmon transcript quantification, IsoformSwitchAnalyzeR import / switch testing, lightweight ISAR visualizations, and multi-contrast summaries. GO enrichment and functional annotation modules are optional because they add extra runtime, storage, and external resource requirements.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

Anton-Bch/isoform-nf-core was originally written by Anton Buch and Tobias Polley.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
