# Anton-Bch/isoform-nf-core: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read QC
- [SRA download](#sra-download) - Optional conversion of SRA runs to FASTQ
- [cat/fastq](#catfastq) - Concatenate multiple sequencing runs per sample
- [fastp](#fastp) - Adapter trimming and read filtering
- [Salmon](#salmon) - Transcriptome indexing and transcript abundance quantification
- [IsoformSwitchAnalyzeR](#isoformswitchanalyzer) - Isoform switch import and differential isoform usage analysis
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### SRA download

<details markdown="1">
<summary>Output files</summary>

- `fetch/`
  - `*.fastq.gz`: FASTQ files downloaded from SRA accessions in `--sra_manifest`.

</details>

When `--sra_manifest` is used instead of `--input`, the pipeline downloads each SRA run with `prefetch`, converts it with `fasterq-dump`, compresses the generated FASTQ files, detects whether each run is single-end or paired-end, and then sends the reads through the same FastQC, concatenation, fastp, Salmon, and IsoformSwitchAnalyzeR path as regular FASTQ input.

### cat/fastq

<details markdown="1">
<summary>Output files</summary>

- `cat/`
  - `*.merged.fastq.gz`: FASTQ files after combining multiple sequencing runs for the same biological sample.

</details>

If the samplesheet contains multiple rows with the same `sample` value, the pipeline first validates that the sample-level metadata are consistent and then concatenates the FASTQ files. This gives downstream tools one read file, or one read pair, per biological sample.

### fastp

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `*.fastp.fastq.gz`: cleaned FASTQ files used for Salmon quantification.
  - `*.fastp.html`: per-sample HTML report.
  - `*.fastp.json`: machine-readable per-sample report used by MultiQC.
  - `*.fastp.log`: fastp command log.

</details>

[fastp](https://github.com/OpenGene/fastp) performs adapter detection, adapter trimming, and quality filtering before transcript quantification. The cleaned reads are the direct input to Salmon.

### Salmon

<details markdown="1">
<summary>Output files</summary>

- `salmon/salmon/`
  - Salmon transcriptome index built from `--transcript_fasta`.
  - If `--genome_fasta` is provided, the index is built as a decoy-aware gentrome index.
- `salmon/<sample>/`
  - `quant.sf`: transcript-level abundance estimates used by IsoformSwitchAnalyzeR.
  - `quant.genes.sf`: gene-level abundance estimates generated using `--gtf`.
  - `cmd_info.json`, `lib_format_counts.json` and `aux_info/`: Salmon run metadata and library-format information.
- `salmon/*_meta_info.json` and `salmon/*_lib_format_counts.json`: selected Salmon metadata files copied to the top-level Salmon output directory for easier MultiQC parsing.

</details>

[Salmon](https://salmon.readthedocs.io/) estimates transcript abundance from the cleaned reads. These transcript-level estimates are the quantitative input for the isoform switch analysis.

### IsoformSwitchAnalyzeR

<details markdown="1">
<summary>Output files</summary>

- `isar/isar_analysis/`
  - `design_matrix.csv`: sample-to-condition table used by IsoformSwitchAnalyzeR.
  - `quant_dirs.csv`: mapping from sample names to Salmon quantification directories.
  - `filtered_annotation.gtf`: GTF subset containing quantified transcripts.
  - `switchAnalyzeRlist_imported.rds`: imported IsoformSwitchAnalyzeR object before statistical testing.
  - `switchAnalyzeRlist_analyzed.rds`: analyzed IsoformSwitchAnalyzeR object after switch testing, when testing is possible.
  - `switch_summary.csv`: summary of detected isoform switches per comparison.
  - `top_switches.csv`: top-ranked switching genes/isoforms.
  - `analysis_notes.txt`: short run summary, including whether the DEXSeq-based switch test was run or skipped.
  - `sessionInfo.txt`: R session information for reproducibility.

</details>

[IsoformSwitchAnalyzeR](https://bioconductor.org/packages/IsoformSwitchAnalyzeR/) imports the Salmon transcript estimates together with the sample design and transcript annotation. If the design contains two conditions with at least two samples per condition, the pipeline runs the DEXSeq-based isoform switch test. Smaller smoke-test datasets can still be imported, but statistical switch testing is skipped and this is recorded in `analysis_notes.txt`.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
