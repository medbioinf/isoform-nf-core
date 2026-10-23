# medbioinf/isoform-nf-core

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/medbioinf/isoform-nf-core)
[![GitHub Actions CI Status](https://github.com/medbioinf/isoform-nf-core/actions/workflows/nf-test.yml/badge.svg)](https://github.com/medbioinf/isoform-nf-core/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/medbioinf/isoform-nf-core/actions/workflows/linting.yml/badge.svg)](https://github.com/medbioinf/isoform-nf-core/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/medbioinf/isoform-nf-core)

## Introduction

**medbioinf/isoform-nf-core** is a reference-based RNA-seq pipeline for isoform switch analysis. It accepts local FASTQ samplesheets, public SRA run manifests, or existing Salmon quantification directories. Starting from reads, it performs QC and preprocessing, quantifies transcript abundance with Salmon, and imports the results into IsoformSwitchAnalyzeR for differential isoform usage analysis. Existing Salmon results enter directly at the IsoformSwitchAnalyzeR stage.

The current workflow runs:

1. Raw read QC with [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
2. Optional SRA download and conversion with [`prefetch` / `fasterq-dump`](https://github.com/ncbi/sra-tools) (SRA Toolkit)
3. Multiple-run concatenation with `cat/fastq`
4. Adapter and quality filtering with [`fastp`](https://github.com/OpenGene/fastp)
5. Transcript quantification with [`Salmon`](https://salmon.readthedocs.io/)
6. Isoform switch import and testing with [`IsoformSwitchAnalyzeR`](https://bioconductor.org/packages/IsoformSwitchAnalyzeR/)
7. ISAR summary visualizations, including per-comparison counts and UpSet-style intersections
8. Optional GO enrichment for significant switching genes with [`WebGestaltR`](https://github.com/bzhanglab/WebGestaltR)
9. Optional functional annotation with [`Pfam`](https://www.ebi.ac.uk/interpro/), [`IUPred2A`](https://iupred2a.elte.hu/), [`SignalP`](https://services.healthtech.dtu.dk/services/SignalP-6.0/), [`DeepTMHMM`](https://dtu.biolib.com/DeepTMHMM), and [`DeepLoc2`](https://services.healthtech.dtu.dk/services/DeepLoc-2.0/)
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
nextflow run medbioinf/isoform-nf-core \
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

medbioinf/isoform-nf-core was originally written by Anton Buch and Tobias Polley.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use this pipeline, please cite the tools it runs. The references below cover the tools used by the pipeline; the same list is also maintained in the [`CITATIONS.md`](CITATIONS.md) file.

### Pipeline framework

- **nf-core**

  > Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

- **Nextflow**

  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. _Nat Biotechnol._ 2017;35(4):316-319. doi: [10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820).

### Read processing and quantification

- **FastQC**

  > Andrews S. FastQC: a quality control tool for high throughput sequence data. 2010. Available online: <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>

- **SRA Toolkit** (`prefetch`, `fasterq-dump`)

  > NCBI SRA Toolkit. National Center for Biotechnology Information. Available online: <https://github.com/ncbi/sra-tools>

- **fastp**

  > Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. _Bioinformatics._ 2018;34(17):i884-i890. doi: [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560).

- **Salmon**

  > Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. Salmon provides fast and bias-aware quantification of transcript expression. _Nat Methods._ 2017;14(4):417-419. doi: [10.1038/nmeth.4197](https://doi.org/10.1038/nmeth.4197).

### Isoform switch analysis and enrichment

- **IsoformSwitchAnalyzeR**

  > Vitting-Seerup K, Sandelin A. IsoformSwitchAnalyzeR: analysis of changes in genome-wide patterns of alternative splicing and its functional consequences. _Bioinformatics._ 2019;35(21):4469-4471. doi: [10.1093/bioinformatics/btz247](https://doi.org/10.1093/bioinformatics/btz247).

  > Vitting-Seerup K, Sandelin A. The landscape of isoform switches in human cancers. _Mol Cancer Res._ 2017;15(9):1206-1220. doi: [10.1158/1541-7786.MCR-16-0459](https://doi.org/10.1158/1541-7786.MCR-16-0459).

- **WebGestaltR**

  > Elizarraras JM, Liao Y, Shi Z, Zhu Q, Pico AR, Zhang B. WebGestalt 2024: faster gene set analysis and new support for metabolomics and multi-omics. _Nucleic Acids Res._ 2024;52(W1):W415-W421. doi: [10.1093/nar/gkae456](https://doi.org/10.1093/nar/gkae456).

### Functional annotation (optional)

- **Pfam**

  > Mistry J, Chuguransky S, Williams L, Qureshi M, Salazar GA, Sonnhammer ELL, Tosatto SCE, Paladin L, Raj S, Richardson LJ, Finn RD, Bateman A. Pfam: the protein families database in 2021. _Nucleic Acids Res._ 2021;49(D1):D412-D419. doi: [10.1093/nar/gkaa913](https://doi.org/10.1093/nar/gkaa913).

- **IUPred2A**

  > Mészáros B, Erdős G, Dosztányi Z. IUPred2A: context-dependent prediction of protein disorder as a function of redox state and protein binding. _Nucleic Acids Res._ 2018;46(W1):W329-W337. doi: [10.1093/nar/gky384](https://doi.org/10.1093/nar/gky384).

- **SignalP 6.0**

  > Teufel F, Almagro Armenteros JJ, Johansen AR, Gíslason MH, Pihl SI, Tsirigos KD, Winther O, Brunak S, von Heijne G, Nielsen H. SignalP 6.0 predicts all five types of signal peptides using protein language models. _Nat Biotechnol._ 2022;40(7):1023-1025. doi: [10.1038/s41587-021-01156-3](https://doi.org/10.1038/s41587-021-01156-3).

- **DeepTMHMM**

  > Hallgren J, Tsirigos KD, Pedersen MD, Almagro Armenteros JJ, Marcatili P, Nielsen H, Krogh A, Winther O. DeepTMHMM predicts alpha and beta transmembrane proteins using deep neural networks. _bioRxiv._ 2022. doi: [10.1101/2022.04.08.487609](https://doi.org/10.1101/2022.04.08.487609).

- **DeepLoc 2.0**

  > Thumuluri V, Almagro Armenteros JJ, Johansen AR, Nielsen H, Winther O. DeepLoc 2.0: multi-label subcellular localization prediction using protein language models. _Nucleic Acids Res._ 2022;50(W1):W228-W234. doi: [10.1093/nar/gkac278](https://doi.org/10.1093/nar/gkac278).

### Reporting

- **MultiQC**

  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. _Bioinformatics._ 2016;32(19):3047-3048. doi: [10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354).

### Software packaging and containerisation

- **Bioconda**

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. _Nat Methods._ 2018;15(7):475-476. doi: [10.1038/s41592-018-0046-7](https://doi.org/10.1038/s41592-018-0046-7).

- **BioContainers**

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. _Bioinformatics._ 2017;33(16):2580-2582. doi: [10.1093/bioinformatics/btx192](https://doi.org/10.1093/bioinformatics/btx192).

- **Docker**

  > Merkel D. Docker: lightweight Linux containers for consistent development and deployment. _Linux Journal._ 2014;2014(239):2.

- **Singularity**

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: scientific containers for mobility of compute. _PLoS ONE._ 2017;12(5):e0177459. doi: [10.1371/journal.pone.0177459](https://doi.org/10.1371/journal.pone.0177459).

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).