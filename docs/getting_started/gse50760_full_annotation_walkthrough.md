# End-To-End Run Guide: GSE50760 With All Annotations

Last updated: 2026-07-17

This guide shows how a researcher can run the current `dev` pipeline from a fresh working directory. It uses the GSE50760 validation subset as a realistic example: 4 primary colorectal cancer samples, 4 matched normal colon samples, and 4 liver metastasis samples. The run performs three pairwise contrasts in one workflow execution and enables all annotation modules except GO enrichment.

The enabled analysis includes:

- Isoform switch testing with IsoformSwitchAnalyzeR
- per-contrast summary plots, including the significant-switch count plot and UpSet-style overlap plot
- Pfam protein domains
- IUPred2A intrinsically disordered regions
- SignalP signal peptides
- DeepTMHMM transmembrane topology
- DeepLoc2 subcellular localization
- annotated switch plots for selected genes
- MultiQC summary report

GO enrichment is intentionally excluded to keep this walkthrough focused on the annotation and visualization path.

## 1. Prepare The Machine

Install these tools on the machine where the pipeline should run:

- Git
- Nextflow
- Docker
- enough disk space for public SRA downloads, Nextflow work directories, containers, and result files

Check the runtime:

```bash
git --version
nextflow -version
docker --version
```

For long runs over SSH, use `tmux` or `screen` so the run continues after disconnects.

### Prepare The Annotation Containers

The Docker profile runs every tool in a container, but not every annotation module handles its container in the same way:

| Module       | What the user must provide                                                                                                                                                                                         |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Pfam         | The pipeline provides the scanning container. The user provides the Pfam database files downloaded in Section 3.                                                                                                   |
| IUPred2A     | The pipeline provides a pinned default image. Because it lacks `ps`, the execution report, timeline, and trace are disabled; provide a compatible `--iupred2a_container` only if those runtime reports are needed. |
| DeepTMHMM    | Nothing beyond a working Docker installation. The pipeline pins its container to an immutable image digest.                                                                                                        |
| SignalP 5.0b | A private or local container containing software obtained under terms suitable for the user's institution.                                                                                                         |
| DeepLoc 2.1  | A private or local container containing software obtained under terms suitable for the user's institution.                                                                                                         |

Docker downloads the pipeline-specified Pfam, IUPred2A, and DeepTMHMM images when they are first needed. Network access to the relevant registries is therefore required unless the images have already been cached. IUPred2A normally uses its default image; see [Prerequisites](prerequisites.md#optional-iupred2a-prerequisites) for the optional override that retains full Nextflow runtime reporting.

SignalP and DeepLoc2 are different because their licenses restrict redistribution. The pipeline cannot provide public default images for them. Before starting this all-annotation walkthrough:

1. Obtain SignalP 5.0b and DeepLoc 2.1 under licenses suitable for the intended use.
2. Build or obtain compatible private container images that expose the expected `signalp` and `deeploc2` commands.
3. Push the images to a private registry accessible from the compute machine, or otherwise make suitable local images available.
4. Run `docker login REGISTRY` first if authentication is required.
5. Replace both placeholder references in the parameter file below with the real image references.

The values `your-private-registry/signalp:5.0b` and `your-private-registry/deeploc:2.1` are placeholders, not downloadable example images. The pipeline deliberately stops during parameter validation when either licensed module is enabled without its corresponding container parameter. See [Prerequisites](prerequisites.md#optional-signalp-prerequisites) for the detailed license and runtime notes.

## 2. Download The Pipeline

Clone the pipeline repository and use the `dev` branch:

```bash
mkdir -p ~/projects
cd ~/projects
git clone https://github.com/Anton-Bch/isoform-nf-core.git
cd isoform-nf-core
git switch dev
git pull --ff-only origin dev
```

If the repository already exists, update it instead:

```bash
cd ~/projects/isoform-nf-core
git fetch origin
git switch dev
git pull --ff-only origin dev
```

## 3. Download Reference Files

Create a local reference directory. The example below uses GENCODE human release 49, matching the validation runs:

```bash
mkdir -p ~/projects/references/gencode_v49
cd ~/projects/references/gencode_v49

wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
```

Download Pfam if Pfam domain annotation should be enabled:

```bash
mkdir -p ~/projects/references/pfam/Pfam37.0
cd ~/projects/references/pfam/Pfam37.0

wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-A.hmm.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-A.hmm.dat.gz
gunzip -k Pfam-A.hmm.gz
gunzip -k Pfam-A.hmm.dat.gz
```

Check that the required files exist:

```bash
test -s ~/projects/references/gencode_v49/gencode.v49.transcripts.fa.gz
test -s ~/projects/references/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz
test -s ~/projects/references/pfam/Pfam37.0/Pfam-A.hmm
test -s ~/projects/references/pfam/Pfam37.0/Pfam-A.hmm.dat
```

## 4. Understand The Metadata Files

The pipeline needs two study-specific metadata files:

- an SRA manifest or FASTQ samplesheet describing the samples
- a contrast file describing which biological groups should be compared

For this GSE50760 example, create both files in a local `metadata/` directory. A researcher running their own experiment would create equivalent files for their own accessions or FASTQs.

### SRA Manifest

Use an SRA manifest when the raw data should be downloaded from public SRA accessions. Create the GSE50760 example manifest:

```bash
cd ~/projects/isoform-nf-core
mkdir -p metadata

cat > metadata/gse50760_4v4v4_sra_manifest.csv <<'EOF'
sample,condition,replicate,run_accession,strandedness,batch
AMC_2-1,primary_crc,1,SRR975551,auto,gse50760
AMC_3-1,primary_crc,2,SRR975552,auto,gse50760
AMC_5-1,primary_crc,3,SRR975553,auto,gse50760
AMC_6-1,primary_crc,4,SRR975554,auto,gse50760
AMC_2-2,normal_colon,1,SRR975569,auto,gse50760
AMC_3-2,normal_colon,2,SRR975570,auto,gse50760
AMC_5-2,normal_colon,3,SRR975571,auto,gse50760
AMC_6-2,normal_colon,4,SRR975572,auto,gse50760
AMC_2-3,liver_metastasis,1,SRR975587,auto,gse50760
AMC_3-3,liver_metastasis,2,SRR975588,auto,gse50760
AMC_5-3,liver_metastasis,3,SRR975589,auto,gse50760
AMC_6-3,liver_metastasis,4,SRR975590,auto,gse50760
EOF

column -s, -t metadata/gse50760_4v4v4_sra_manifest.csv | sed -n '1,15p'
```

The important idea is one row per sequencing run, with a stable sample name and the biological condition/group for that sample. For this example, the condition values are:

```text
normal_colon
primary_crc
liver_metastasis
```

Those condition names are later reused in the contrast file and in the annotated switch plot parameters.

The example SRA manifest has this structure:

```csv
sample,condition,replicate,run_accession,strandedness,batch
AMC_2-1,primary_crc,1,SRR975551,auto,gse50760
AMC_3-1,primary_crc,2,SRR975552,auto,gse50760
```

For a new public dataset, replace the `sample`, `condition`, `replicate`, and `run_accession` values with the new study metadata. `strandedness` can be set to `auto` when the library type should be inferred. `batch` is optional study information that helps keep metadata readable.

### Contrast File

The contrast file tells the pipeline which group comparisons to test in IsoformSwitchAnalyzeR. Create the GSE50760 example contrast file:

```bash
cat > metadata/gse50760_4v4v4_multicontrasts.csv <<'EOF'
contrast,case,control
primary_crc_vs_normal_colon,primary_crc,normal_colon
liver_metastasis_vs_normal_colon,liver_metastasis,normal_colon
liver_metastasis_vs_primary_crc,liver_metastasis,primary_crc
EOF

column -s, -t metadata/gse50760_4v4v4_multicontrasts.csv
```

The GSE50760 validation subset runs three comparisons in one pipeline run:

```csv
contrast,case,control
primary_crc_vs_normal_colon,primary_crc,normal_colon
liver_metastasis_vs_normal_colon,liver_metastasis,normal_colon
liver_metastasis_vs_primary_crc,liver_metastasis,primary_crc
```

For another study, create one row per comparison. The condition names must match the condition/group column in the sample metadata. This is how a study with more than two groups can be analyzed without launching separate pipeline runs for every pairwise comparison.

## 5. Create The Parameter File

Return to the pipeline repository and create a dated parameter file. Keeping the date at the beginning of the output name makes repeated validation runs easier to sort:

```bash
cd ~/projects/isoform-nf-core
mkdir -p params

cat > params/20260705_gse50760_4v4v4_all_annotations.yml <<EOF
sra_manifest: "metadata/gse50760_4v4v4_sra_manifest.csv"
contrasts: "metadata/gse50760_4v4v4_multicontrasts.csv"

transcript_fasta: "${HOME}/projects/references/gencode_v49/gencode.v49.transcripts.fa.gz"
gtf: "${HOME}/projects/references/gencode_v49/gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
pfam_db: "${HOME}/projects/references/pfam/Pfam37.0"

run_isar: true
run_isar_visualization: true
run_isar_contrast_summary: true

run_iupred2a: true
run_signalp: true
signalp_container: "your-private-registry/signalp:5.0b"
run_deeptmhmm: true
run_deeploc2: true
deeploc2_container: "your-private-registry/deeploc:2.1"

run_annotated_switch_plots: true
annotated_switch_genes: "ZNRF3,PBX3,YEATS4"
annotated_switch_condition1: "normal_colon"
annotated_switch_condition2: "primary_crc"

outdir: "results/20260705_gse50760_4v4v4_all_annotations"
multiqc_title: "20260705_gse50760_4v4v4_all_annotations"
EOF
```

Adjust the absolute reference paths if the files were downloaded somewhere else. The annotated switch plot settings above request gene-level plots for `ZNRF3`, `PBX3`, and `YEATS4` in the normal-colon versus primary-CRC comparison. The full ISAR analysis still evaluates all three contrasts from the contrast file.

Replace `your-private-registry/signalp:5.0b` and `your-private-registry/deeploc:2.1` with the licensed images prepared above before launching the workflow. Do not run the example with the placeholder values unchanged.

## 6. Start The Run

Start a detachable terminal session:

```bash
tmux new -s gse50760_4v4v4_all_annotations
```

Inside `tmux`, launch the workflow:

```bash
cd ~/projects/isoform-nf-core

nextflow run . \
    -profile docker \
    -params-file params/20260705_gse50760_4v4v4_all_annotations.yml \
    -work-dir work_20260705_gse50760_4v4v4_all_annotations
```

Detach from `tmux` with `Ctrl-b`, then `d`. Reattach later with:

```bash
tmux attach -t gse50760_4v4v4_all_annotations
```

If the run stops because of a temporary problem, fix the issue and restart with `-resume`:

```bash
nextflow run . \
    -profile docker \
    -params-file params/20260705_gse50760_4v4v4_all_annotations.yml \
    -work-dir work_20260705_gse50760_4v4v4_all_annotations \
    -resume
```

## 7. Monitor Progress

From the pipeline repository:

```bash
nextflow log
tail -f .nextflow.log
```

Check running containers and disk use:

```bash
docker ps
df -h .
du -sh work_20260705_gse50760_4v4v4_all_annotations results/20260705_gse50760_4v4v4_all_annotations 2>/dev/null
```

## 8. Understand The Output Structure

The final results are written to:

```text
results/20260705_gse50760_4v4v4_all_annotations/
```

The most useful starting point is:

```text
multiqc/multiqc_report.html
```

Important result folders:

```text
fastqc/                         raw-read quality control
fastp/                          read trimming/filtering reports
salmon/                         transcript quantification
isar/isar_analysis/             IsoformSwitchAnalyzeR objects and tables
isar/isar_visualization/        one general switch-analysis folder per contrast
isar/isar_contrast_summary/     multi-contrast summary plots and tables
pfam/                           protein domain annotation and summaries
iupred2a/                       intrinsically disordered region annotation
signalp/                        signal peptide annotation
deeptmhmm/                      transmembrane topology annotation
deeploc2/                       subcellular localization annotation
annotated/annotated_switch_plots/ selected gene-level annotated switch plots
multiqc/                        combined HTML report
pipeline_info/                  execution reports, trace files, and software versions
```

Key files to inspect first:

```text
isar/isar_contrast_summary/significant_isoform_switches_per_comparison.csv
isar/isar_contrast_summary/significant_isoform_switches_per_comparison.png
isar/isar_contrast_summary/isoform_switch_upset.png
isar/isar_visualization/comparison_visualizations.csv
isar/isar_visualization/primary_crc_vs_normal_colon/isoform_switch_volcano.png
isar/isar_visualization/primary_crc_vs_normal_colon/top_switching_genes.png
annotated/annotated_switch_plots/annotation_status.csv
annotated/annotated_switch_plots/annotated_switch_plot_summary.csv
annotated/annotated_switch_plots/*_annotated_switch.png
multiqc/multiqc_report.html
```

The `annotation_status.csv` file is the quickest check for whether Pfam, IUPred2A, SignalP, DeepTMHMM, and DeepLoc2 were available to the annotated switch plot step.

## 9. Quick Success Checks

After the workflow finishes:

```bash
OUT=results/20260705_gse50760_4v4v4_all_annotations

test -s "$OUT/isar/isar_analysis/switchAnalyzeRlist_analyzed.rds"
test -s "$OUT/isar/isar_contrast_summary/isoform_switch_upset.png"
test -s "$OUT/isar/isar_visualization/primary_crc_vs_normal_colon/isoform_switch_volcano.png"
test -s "$OUT/pfam/pfam_import/pfam_domain_summary.csv"
test -s "$OUT/iupred2a/iupred2a_import/iupred2a_idr_summary.csv"
test -s "$OUT/signalp/signalp_import/signalp_summary.csv"
test -s "$OUT/deeptmhmm/deeptmhmm_import/deeptmhmm_summary.csv"
test -s "$OUT/deeploc2/deeploc2_import/deeploc2_summary.csv"
test -s "$OUT/annotated/annotated_switch_plots/annotation_status.csv"

column -s, -t "$OUT/annotated/annotated_switch_plots/annotation_status.csv"
```

If one annotation layer is unavailable, inspect the matching module folder and the files in `pipeline_info/` before rerunning. In many cases, `-resume` will reuse all completed upstream work.

## 10. What To Copy Back

For sharing or manual inspection, copy the final `results/20260705_gse50760_4v4v4_all_annotations/` directory. If bandwidth or storage is limited, copy these first:

```text
results/20260705_gse50760_4v4v4_all_annotations/multiqc/
results/20260705_gse50760_4v4v4_all_annotations/isar/
results/20260705_gse50760_4v4v4_all_annotations/annotated/
results/20260705_gse50760_4v4v4_all_annotations/pfam/
results/20260705_gse50760_4v4v4_all_annotations/iupred2a/
results/20260705_gse50760_4v4v4_all_annotations/signalp/
results/20260705_gse50760_4v4v4_all_annotations/deeptmhmm/
results/20260705_gse50760_4v4v4_all_annotations/deeploc2/
results/20260705_gse50760_4v4v4_all_annotations/pipeline_info/
```

The `work_20260705_gse50760_4v4v4_all_annotations/` directory should usually not be copied back. It is large and mainly useful for `-resume` on the same machine.

Example copy command from a local machine:

```bash
rsync -avP ubuntu@YOUR_VM:/home/ubuntu/projects/isoform-nf-core/results/20260705_gse50760_4v4v4_all_annotations/ \
    ./20260705_gse50760_4v4v4_all_annotations/
```

After the results have been inspected and no resume is needed, the work directory can be removed on the compute machine:

```bash
rm -rf work_20260705_gse50760_4v4v4_all_annotations
```
