# Anton-Bch/isoform-nf-core: Documentation

Last updated: 2026-07-05

This documentation is grouped by use case. Start with the first section if you are new to Nextflow, nf-core, or isoform switch analysis.

## Start Here

- [Project structure and nf-core basics](getting_started/project_structure_and_nfcore_basics.md): repository layout, Nextflow concepts, nf-core conventions, modules, configuration, and outputs.
- [Tools and biology explained](getting_started/tools_and_biology_explained.md): non-biologist explanation of genes, transcripts, isoforms, read types, references, and the tools used in the workflow.
- [Interface and inputs explained](getting_started/interface_and_inputs_explained.md): FASTQ samplesheets, SRA manifests, contrasts, reference files, user parameters, and validation rules.
- [Prerequisites to run the pipeline](getting_started/prerequisites.md): runtime setup, input/reference requirements, optional Pfam resources, and compute/storage expectations.

## User Docs

- [Usage](usage.md): how to run the pipeline and provide inputs.
- [Output](output.md): result folders, key files, and how to interpret outputs.

## Developer Docs

- [Workflow implementation walkthrough](developer/workflow_implementation_walkthrough.md): code-oriented walkthrough of `main.nf`, `workflows/isoform.nf`, modules, channels, helper scripts, and extension points.

## Operations

- [Runtime monitoring](operations/runtime_monitoring.md): lightweight VM monitoring notes for long-running Nextflow executions.
- [Implementation validation log](operations/implementation_validation_log.md): concrete checks performed for FASTQ, SRA, ISAR, visualization, annotation, and multi-contrast implementation steps.
