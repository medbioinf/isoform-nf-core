# Anton-Bch/isoform-nf-core: Output

## Introduction

This document describes the output produced by the pipeline. MultiQC summarizes run-level QC, while the ISAR, contrast-summary, annotation, and annotated-switch-plot modules write their own tables and plots into dedicated result directories.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [FastQC](#fastqc) - Raw read QC
- [SRA download](#sra-download) - Optional conversion of SRA runs to FASTQ
- [cat/fastq](#catfastq) - Concatenate multiple sequencing runs per sample
- [fastp](#fastp) - Adapter trimming and read filtering
- [Salmon](#salmon) - Transcriptome indexing and transcript abundance quantification
- [IsoformSwitchAnalyzeR](#isoformswitchanalyzer) - Isoform switch import and differential isoform usage analysis
- [ISAR visualization](#isar-visualization) - Lightweight plots and candidate tables from IsoformSwitchAnalyzeR output
- [ISAR contrast summary](#isar-contrast-summary) - Cross-comparison bar and UpSet-style plots for significant isoform switches
- [GO enrichment](#go-enrichment) - Optional gene ontology enrichment for significant switching genes
- [Pfam annotation](#pfam-annotation) - Optional protein-domain annotation for significant switch candidates
- [IUPred2A annotation](#iupred2a-annotation) - Optional intrinsically disordered region and ANCHOR2 annotation
- [SignalP annotation](#signalp-annotation) - Optional signal peptide annotation
- [DeepTMHMM annotation](#deeptmhmm-annotation) - Optional transmembrane topology annotation
- [DeepLoc2 annotation](#deeploc2-annotation) - Optional subcellular localization annotation
- [Annotated switch plots](#annotated-switch-plots) - Optional IsoformSwitchAnalyzeR switch plots with available annotation tracks
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
- `salmon/*_requested_lib_type.txt`: the library type requested by the pipeline (`A` records true Salmon automatic inference); the inferred type remains available in `*_lib_format_counts.json`.

</details>

[Salmon](https://salmon.readthedocs.io/) estimates transcript abundance from the cleaned reads. These transcript-level estimates are the quantitative input for the isoform switch analysis. When `--salmon_input` is used, this entire pipeline stage is skipped and no new `salmon/` output is published; ISAR consumes the supplied quantification directories instead.

### IsoformSwitchAnalyzeR

<details markdown="1">
<summary>Output files</summary>

- `isar/isar_analysis/`
  - `design_matrix.csv`: sample-to-condition table used by IsoformSwitchAnalyzeR.
  - `quant_dirs.csv`: mapping from sample names to Salmon quantification directories.
  - `comparisons.csv`: requested pairwise comparisons and their labels, when statistical testing is configured.
  - `filtered_annotation.gtf`: GTF subset containing quantified transcripts.
  - `switchAnalyzeRlist_imported.rds`: imported IsoformSwitchAnalyzeR object before statistical testing.
  - `switchAnalyzeRlist_analyzed.rds`: analyzed IsoformSwitchAnalyzeR object after switch testing, when testing is possible.
  - `switch_summary.csv`: summary of detected isoform switches per comparison.
  - `top_switches.csv`: top-ranked switching genes/isoforms.
  - `go_gene_scores.csv`: gene-level GO-enrichment handoff table with per-contrast q-value, dIF, and significance columns.
  - `analysis_notes.txt`: short run summary, including whether the DEXSeq-based switch test was run or skipped.
  - `sessionInfo.txt`: R session information for reproducibility.

</details>

[IsoformSwitchAnalyzeR](https://bioconductor.org/packages/IsoformSwitchAnalyzeR/) imports the Salmon transcript estimates together with the sample design and transcript annotation. If a contrast file is supplied, the pipeline runs the requested pairwise comparisons in one ISAR analysis when each contrast condition has at least two samples. If no contrast file is supplied, a two-condition dataset is tested as a single comparison. Smaller smoke-test datasets can still be imported, but statistical switch testing is skipped and this is recorded in `analysis_notes.txt`.

### ISAR visualization

<details markdown="1">
<summary>Output files</summary>

- `isar/isar_visualization/`
  - `comparison_visualizations.csv`: comparison names, condition pairs, output directories, and completion status.
  - `visualization_notes.txt`: run-level summary and explanation if plots were skipped.
  - `<comparison>/top_switch_plots.pdf`: official IsoformSwitchAnalyzeR `switchPlot()` pages for that comparison's top-ranked switching genes.
  - `<comparison>/isoform_switch_volcano.png` / `.pdf`: isoform-level effect size versus statistical support for one comparison. Labels identify the top-N significant isoform points by q-value; repeated gene labels indicate that multiple isoforms from that gene are among the top points.
  - `<comparison>/top_switching_genes.png` / `.pdf`: top genes ranked within that comparison by gene-level switch q-value.
  - `<comparison>/top_gene_isoform_usage.pdf`: combined PDF with per-gene isoform usage plots for that comparison.
  - `<comparison>/top_gene_isoform_usage/*.png` / `.pdf`: individual per-gene isoform usage plots.
  - `<comparison>/ptc_switch_summary.png` / `.pdf` / `.csv`: PTC status summary for significant switch candidates.
  - `<comparison>/top_isoform_candidates.csv`: top isoform-level candidates from that comparison.
  - `<comparison>/top_gene_summary.csv`: top gene-level summary for that comparison.
  - `<comparison>/visualization_notes.txt`: comparison-specific run summary.

</details>

The ISAR visualization step is controlled by `--run_isar_visualization` and `--isar_visualization_top_n`. Every tested comparison receives a separate subdirectory, volcano plot, and independently ranked top-N set. It is intentionally lightweight and does not require external annotation databases. If the ISAR step only imported data, or if no tested comparisons are available, the visualization step exits successfully and writes explanatory notes instead of failing the workflow.

### ISAR contrast summary

<details markdown="1">
<summary>Output files</summary>

- `isar/isar_contrast_summary/`
  - `significant_isoform_switches_per_comparison.csv`: significant isoform-switch counts per comparison.
  - `significant_isoform_switches_per_comparison.png` / `.pdf`: bar plot of significant isoform switches per comparison.
  - `isoform_switch_intersections.csv`: largest UpSet-style isoform-switch intersections.
  - `isoform_switch_intersection_members.csv`: isoforms and genes contributing to each plotted intersection.
  - `isoform_switch_upset.png` / `.pdf`: UpSet-style intersection plot with left-side set-size bars for significant isoforms per comparison.
  - `isar_contrast_summary_notes.txt`: summary of inputs and any skipped-output reasons.

</details>

The contrast summary is generated from the analyzed IsoformSwitchAnalyzeR result and uses the configured `--isar_qvalue_cutoff` and `--isar_dif_cutoff`. It is most informative when `--contrasts` contains multiple pairwise comparisons.

### GO enrichment

GO enrichment support is optional and controlled by `--run_go_enrichment`. It summarizes genes with significant isoform switches per contrast and runs WebGestaltR over-representation analysis against Gene Ontology or another configured WebGestaltR database.

- `go/`
  - `go_input/go_gene_scores.csv`: copy of the ISAR-generated per-contrast gene table used as GO enrichment input.
  - `go_input/*_significant_genes.txt`: significant switching genes submitted for each contrast.
  - `go_input/*_background_genes.txt`: background genes used for each contrast.
  - `go_enrichment.csv`: combined enrichment table across contrasts.
  - `*_go_enrichment.csv`: per-contrast enrichment result tables.
  - `go_summary.txt`: status, gene counts, database settings, and any skip/failure notes.
  - `webgestalt_report/`: WebGestaltR report directories when enrichment is executed.

The GO module consumes `isar/isar_analysis/go_gene_scores.csv` rather than the serialized ISAR R object, so it only needs WebGestaltR at runtime. The default background universe is all genes tested by IsoformSwitchAnalyzeR for the corresponding contrast. Use `--go_reference_gene_file` only when a custom newline-delimited background gene universe is needed. GO enrichment is gene-level, so it complements rather than replaces isoform-level annotated switch plots.

### Pfam annotation

Pfam support is optional and split into three parts. First, the pipeline prepares amino-acid FASTA input from significant ISAR switch candidates. Second, it either runs `pfam_scan.pl` itself when `--pfam_db` is provided, or imports an existing `pfam_scan.pl` result file when `--pfam_results` is provided. Third, it imports those domain hits into IsoformSwitchAnalyzeR and creates domain consequence plots.

- `pfam/pfam_prepare/`
  - `isoform_pfam_candidates_AA.fasta`: amino-acid FASTA for external `pfam_scan.pl` / HMMER scanning.
  - `pfam_top_candidates.csv`: strongest isoform switch candidates selected for Pfam inspection.
  - `switchAnalyzeRlist_with_sequences.rds`: ISAR object after sequence extraction, when sequence extraction succeeds.
  - `pfam_prepare_notes.txt`: summary of candidate and FASTA generation.
- `pfam/pfam_scan.out`: `pfam_scan.pl` output generated by the pipeline when `--pfam_db` is provided and `--pfam_results` is not provided.
- `pfam/pfam_scan.log`: Pfam scan run log.
- `pfam/pfam_import/`
  - `switchAnalyzeRlist_with_pfam.rds`: ISAR object after `analyzePFAM()`.
  - `pfam_domain_analysis.csv`: imported Pfam domain hit table.
  - `isoform_features_with_pfam.csv`: isoform feature table with Pfam annotation columns.
  - `significant_switches_with_pfam.csv`: significant switch candidates after Pfam import.
  - `pfam_switch_consequences.csv`: domain gain/loss/switch consequences from IsoformSwitchAnalyzeR.
  - `pfam_domain_summary.png` / `.pdf` / `.csv`: overview of domain annotation availability.
  - `pfam_import_notes.txt`: import summary and warnings.
- `pfam/pfam_visualization/`
  - `pfam_switch_consequence_summary.png` / `.pdf` / `.csv`: counts of domain consequences, grouped by usage-shift direction.
  - `pfam_top_domain_architectures.pdf`: combined PDF of top domain architecture changes.
  - `domain_architecture/*.png` / `.pdf`: individual domain architecture plots.
  - `pfam_top_domain_change_candidates.csv`: ranked domain-change candidates.
  - `pfam_visualization_notes.txt`: interpretation notes.

The pipeline does not download the Pfam database automatically. The Pfam database is large, so it should be prepared as a reference resource and passed with `--pfam_db`.

### IUPred2A annotation

IUPred2A support is optional and controlled by `--run_iupred2a`. It predicts intrinsically disordered regions (IDRs) and ANCHOR2 binding regions in proteins from significant isoform switch candidates.

- `iupred2a/iupred2a_prepare/`
  - `isoform_iupred2a_candidates_AA.fasta`: amino-acid FASTA used as IUPred2A input.
  - `iupred2a_top_candidates.csv`: strongest isoform switch candidates selected for IUPred2A inspection.
  - `switchAnalyzeRlist_with_sequences.rds`: ISAR object after sequence extraction, when sequence extraction succeeds.
  - `iupred2a_prepare_notes.txt`: summary of candidate and FASTA generation.
- `iupred2a/iupred2a/`
  - `iupred2a_anchor2.raw`: raw IUPred2A/ANCHOR2 output.
  - `iupred2a_anchor2_isar.out`: converted output in the block format imported by IsoformSwitchAnalyzeR.
  - `iupred2a.log`: run and conversion log.
- `iupred2a/iupred2a_import/`
  - `switchAnalyzeRlist_with_iupred2a.rds`: ISAR object after `analyzeIUPred2A()`.
  - `iupred2a_idr_analysis.csv`: imported IDR annotation table.
  - `isoform_features_with_iupred2a.csv`: isoform feature table with IDR annotation columns when available.
  - `iupred2a_idr_summary.png` / `.pdf` / `.csv`: overview of imported IDR annotation availability.
  - `iupred2a_import_notes.txt`: import summary and interpretation notes.

### SignalP annotation

SignalP support is optional and controlled by `--run_signalp`. It predicts N-terminal signal peptides in proteins from significant isoform switch candidates. A positive prediction supports entry into the secretory pathway, but does not by itself establish final secretion, localization, or a membrane anchor. A licensed user-supplied SignalP 5.0b container is required via `--signalp_container`.

- `signalp/signalp_prepare/`
  - `isoform_signalp_candidates_AA.fasta`: amino-acid FASTA used as SignalP input.
  - `signalp_top_candidates.csv`: strongest isoform switch candidates selected for SignalP inspection.
  - `switchAnalyzeRlist_with_sequences.rds`: ISAR object after sequence extraction, when sequence extraction succeeds.
  - `signalp_prepare_notes.txt`: summary of candidate and FASTA generation.
- `signalp/signalp/`
  - `signalp5_summary.signalp5`: SignalP 5 short-format output imported by IsoformSwitchAnalyzeR.
  - `signalp.log`: SignalP run log.
- `signalp/signalp_import/`
  - `switchAnalyzeRlist_with_signalp.rds`: ISAR object after `analyzeSignalP()`.
  - `signalp_signal_peptide_analysis.csv`: imported signal peptide annotation table.
  - `isoform_features_with_signalp.csv`: isoform feature table with signal peptide annotation columns when available.
  - `signalp_summary.png` / `.pdf` / `.csv`: overview of imported signal peptide annotation availability.
  - `signalp_import_notes.txt`: import summary and interpretation notes.

The pipeline does not redistribute or automatically download SignalP. Use a container provisioned under the appropriate institutional license and record an immutable image digest or SIF checksum for reproducibility.

### DeepTMHMM annotation

DeepTMHMM support is optional and controlled by `--run_deeptmhmm`. It predicts transmembrane topology regions in proteins from significant isoform switch candidates.

- `deeptmhmm/deeptmhmm_prepare/`
  - `isoform_deeptmhmm_candidates_AA.fasta`: amino-acid FASTA used as DeepTMHMM input.
  - `deeptmhmm_top_candidates.csv`: strongest isoform switch candidates selected for inspection.
  - `switchAnalyzeRlist_with_sequences.rds`: ISAR object after sequence extraction, when sequence extraction succeeds.
  - `deeptmhmm_prepare_notes.txt`: summary of candidate and FASTA generation.
- `deeptmhmm/deeptmhmm/`
  - `predicted_topologies.3line`: raw DeepTMHMM three-line topology output.
  - `deeptmhmm_regions_isar.tsv`: converted topology regions imported by IsoformSwitchAnalyzeR.
  - `deeptmhmm.log`: DeepTMHMM run log.
- `deeptmhmm/deeptmhmm_import/`
  - `switchAnalyzeRlist_with_deeptmhmm.rds`: ISAR object after `analyzeDeepTMHMM()`.
  - `deeptmhmm_topology_analysis.csv`: imported topology region table.
  - `isoform_features_with_deeptmhmm.csv`: isoform feature table after topology import.
  - `deeptmhmm_summary.png` / `.pdf` / `.csv`: overview of topology annotation availability.
  - `deeptmhmm_import_notes.txt`: import summary and interpretation notes.

DeepTMHMM can be slow on real datasets because it embeds protein sequences before topology prediction. Keep it optional for routine quantification-only runs.

### DeepLoc2 annotation

DeepLoc2 support is optional and controlled by `--run_deeploc2`. It predicts subcellular localization labels for proteins from significant isoform switch candidates. A licensed user-supplied DeepLoc 2.1 container is required via `--deeploc2_container`.

- `deeploc2/deeploc2_prepare/`
  - `isoform_deeploc2_candidates_AA.fasta`: amino-acid FASTA used as DeepLoc2 input.
  - `deeploc2_top_candidates.csv`: strongest isoform switch candidates selected for inspection.
  - `switchAnalyzeRlist_with_sequences.rds`: ISAR object after sequence extraction, when sequence extraction succeeds.
  - `deeploc2_prepare_notes.txt`: summary of candidate and FASTA generation.
- `deeploc2/deeploc2/`
  - `results_*.csv`: raw DeepLoc2 prediction table.
  - `deeploc2_isar.csv`: converted prediction table imported by IsoformSwitchAnalyzeR.
  - `deeploc2.log`: DeepLoc2 run log.
- `deeploc2/deeploc2_import/`
  - `switchAnalyzeRlist_with_deeploc2.rds`: ISAR object after `analyzeDeepLoc2()`.
  - `deeploc2_location_analysis.csv`: imported localization table.
  - `isoform_features_with_deeploc2.csv`: isoform feature table after localization import.
  - `deeploc2_summary.png` / `.pdf` / `.csv`: overview of imported localization labels.
  - `deeploc2_import_notes.txt`: import summary and interpretation notes.

The module uses the fast CPU model. Compared with the slower model this trades a small amount of accuracy for throughput, and DeepLoc truncates sequences longer than 1022 amino acids by retaining their ends. Preserve the raw prediction table and run log when reporting results.

### Annotated switch plots

Annotated switch plots are optional and controlled by `--run_annotated_switch_plots`. This step uses the newest available annotated IsoformSwitchAnalyzeR object from the optional annotation chain and renders `switchPlot()` pages for selected switching genes.

- `annotated/annotated_switch_plots/`
  - `annotated_switch_plot_comparisons.csv`: comparison names, output directories, and completion status when all comparisons are plotted.
  - `annotated_switch_plot_notes.txt`: run-level summary.
  - `<comparison>/*_annotated_switch.png` / `.pdf`: per-gene annotated switch plots for one comparison.
  - `<comparison>/annotated_switch_plot_genes.csv`: significant switching genes selected independently for that comparison.
  - `<comparison>/annotated_switch_plot_summary.csv`: plotted genes, strongest isoforms, q-values, dIF values, and output filenames.
  - `<comparison>/annotation_status.csv`: annotation layers available to that comparison's plots.
  - `<comparison>/annotated_switch_plot_notes.txt`: comparison-specific summary and interpretation notes.

By default, every tested comparison receives its own directory and its own top-N ranking. Setting both `--annotated_switch_condition1` and `--annotated_switch_condition2` restricts output to that comparison; in this explicit single-comparison mode, the files are written directly in `annotated_switch_plots/`.

The available tracks depend on which optional annotation modules were run before this step. ORF/PTC information comes from the base ISAR object. Pfam, SignalP, IUPred2A, DeepTMHMM, and DeepLoc2 tracks appear only when their corresponding annotations have been imported.

If `--annotated_switch_genes` is supplied, enabled annotation preparation modules force-include those genes in their protein FASTA target sets. This makes requested gene-level plots more reproducible across runs where the automatic significant-switch candidate list differs slightly.

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
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`. The execution report, timeline, and trace are omitted when the default IUPred2A container is used because that image lacks the process-monitoring utility `ps`; the DAG remains available.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
