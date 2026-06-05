# Anton-Bch/isoform-nf-core: Documentation

Last updated: 2026-06-05

This documentation is grouped by use case. Start with the first section if you are new to Nextflow, nf-core, or isoform switch analysis.

## Start Here

- [Project structure and nf-core basics](getting_started/project_structure_and_nfcore_basics.md): repository layout, Nextflow concepts, nf-core conventions, modules, configuration, and outputs.
- [Tools and biology explained](getting_started/tools_and_biology_explained.md): non-biologist explanation of genes, transcripts, isoforms, read types, references, and the tools used in the workflow.
- [Interface and inputs explained](getting_started/interface_and_inputs_explained.md): FASTQ samplesheets, SRA manifests, contrasts, reference files, user parameters, and validation rules.

## User Docs

- [Usage](usage.md): how to run the pipeline and provide inputs.
- [Output](output.md): result folders, key files, and how to interpret outputs.

## Developer Docs

- [Workflow implementation walkthrough](developer/workflow_implementation_walkthrough.md): code-oriented walkthrough of `main.nf`, `workflows/isoform.nf`, modules, channels, helper scripts, and extension points.

## Design and Roadmap

- [V1 interface status](design/v1_interface_design.md): current implemented input contract and important V1 limitations.
- [V1 architecture status](design/v1_architecture_plan.md): current module structure and remaining architecture decisions.
- [V1 testing status and strategy](design/v1_testing_strategy.md): synthetic and real-data validation status plus remaining test work.
- [V1 packaging and roadmap status](design/v1_packaging_and_roadmap.md): current runtime/container choices and remaining hardening tasks.

## Operations

- [Runtime monitoring](operations/runtime_monitoring.md): lightweight VM monitoring notes for long-running Nextflow executions.
