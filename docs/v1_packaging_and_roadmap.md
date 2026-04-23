# V1 Packaging And Roadmap

This document combines the remaining high-level implementation guidance for V1:

- software packaging and container strategy
- the milestone roadmap for moving from scaffold to working pipeline

## Packaging strategy

The packaging approach should follow a simple rule:

- use shared package/container definitions for standard tools
- keep custom packaging only where the pipeline is genuinely unique

## Standard tools

For common workflow steps, prefer shared `nf-core` modules with their normal package and container definitions.

This applies to:

- `fastqc`
- `fastp`
- `multiqc`
- `salmon/index`
- `salmon/quant`

Why:

- lower maintenance burden
- better alignment with `nf-core`
- easier linting and testing
- less custom Docker maintenance

## ISAR packaging

The preferred first implementation for the `IsoformSwitchAnalyzeR` step is:

- a local module
- package-based environment definitions
- no custom Dockerfile at first

Why this is realistic:

- `IsoformSwitchAnalyzeR` is available through Bioconda
- important dependencies such as `DEXSeq`, `tximport`, and `pfamAnalyzeR` are also packaged
- the current V1 scope does not require every optional external annotation tool in the ISAR vignette

This means the first implementation should try:

- `bioconductor-isoformswitchanalyzer`
- `r-optparse`

plus any additional package dependencies required by the local script

## When to fall back to a custom container

A dedicated custom container for the ISAR step is still acceptable if:

- dependency solving becomes unstable
- package-based execution is too fragile across environments
- extra external annotation tools are added later and become difficult to package cleanly

If that happens, keep the custom container limited to:

- the ISAR module only

Do not use one large pipeline-wide custom image.

## SRA packaging

The optional SRA mode can be implemented in two ways:

1. reuse shared downloader-related `nf-core` components if they fit well
2. keep local SRA download / conversion modules if that is simpler and more reliable

For V1, correctness matters more than maximal reuse.

If local SRA modules produce a cleaner workflow, that is acceptable.

## What should stay out of V1 packaging

Do not try to package every optional ISAR companion tool into the first release.

Examples:

- `CPAT`
- `CPC2`
- `Pfam`
- `SignalP`
- `NetSurfP-2`
- `IUPred2A`
- `DeepLoc2`
- `DeepTMHMM`

These should only be considered later if:

- they answer a real biological question
- they can be run reproducibly offline
- they do not make the pipeline fragile

This includes requests to integrate tools such as `Pfam` immediately just because they are mentioned in the `IsoformSwitchAnalyzeR` vignette.

For V1, the default plan should remain:

- implement the core ISAR path first
- add advanced companion tools later only if clearly justified

## Competitor / landscape analysis

Alongside implementation, the project should keep a lightweight competitor analysis in mind.

Purpose:

- understand how similar pipelines present outputs
- understand which workflow pieces are already standard
- understand what makes this pipeline distinct

Relevant comparison points include:

- `nf-core/rnaseq`
- transcript-level quantification workflows built around `Salmon`
- existing `IsoformSwitchAnalyzeR`-based workflows or wrappers

This comparison should inform positioning and documentation, but it should not block V1 implementation.

## Roadmap overview

The roadmap should move in small, testable milestones.

## Milestone 1: schema and parameter baseline

Goal:

- replace template placeholders with the real V1 interface

Tasks:

- update `nextflow.config`
- update `nextflow_schema.json`
- replace the template input schema with the V1 FASTQ schema
- add `sra_manifest` to the schema
- define the first real user-facing parameters

Done when:

- the new schema matches `docs/v1_interface_design.md`
- help output reflects the real pipeline contract

## Milestone 2: install shared modules

Goal:

- bring in the standard building blocks

Tasks:

- install `fastp`
- install `salmon/index`
- install `salmon/quant`

Done when:

- the repo contains the required shared modules
- `conf/modules.config` is updated as needed

## Milestone 3: input handling

Goal:

- support both input modes with one normalized internal contract

Tasks:

- implement FASTQ input parsing
- implement SRA manifest parsing
- normalize to one internal sample representation
- emit one metadata CSV / channel shape for downstream analysis

Done when:

- both modes feed the same downstream subworkflows

## Milestone 4: preprocessing and quantification

Goal:

- make the generic RNA-seq path work end to end without ISAR yet

Tasks:

- wire `FASTQC`
- wire optional `FASTP`
- wire `SALMON_INDEX`
- wire `SALMON_QUANT`
- wire `MULTIQC`

Done when:

- paired-end and single-end smoke tests can run through Salmon successfully

## Milestone 5: local ISAR module

Goal:

- make the core biology-specific analysis step work in the new repo

Tasks:

- port the reusable R analysis logic
- create a local `ISAR_ANALYSIS` module
- define module inputs and outputs cleanly
- connect the module into the workflow

Done when:

- minimal ISAR outputs are produced from synthetic test runs
- the pipeline’s scope is clearly documented as reference-based isoform-switch analysis rather than novel isoform discovery

## Milestone 6: tests

Goal:

- replace all template test placeholders

Tasks:

- build paired-end synthetic smoke test
- build single-end synthetic smoke test
- update `conf/test.config`
- add `conf/test_single.config`
- update `tests/default.nf.test`
- add an ISAR-enabled minimal test later if stable

Done when:

- the pipeline has meaningful `nf-test` coverage for both read layouts

## Milestone 7: docs and citations

Goal:

- make the repo understandable and publishable

Tasks:

- replace template README text
- document the V1 input formats
- document outputs
- update `CITATIONS.md`
- acknowledge reused `nf-core` components
- make the main outputs easy to understand and present in project discussions
- summarize the pipeline’s scope and limitations clearly, including the lack of de novo novel isoform discovery in V1

Done when:

- the repo docs describe the real pipeline, not the template

## Milestone 8: lint and refinement

Goal:

- bring the implementation into good `nf-core` shape

Tasks:

- run `nf-core pipelines lint`
- fix schema, docs, config, and component issues
- refine naming, labels, and resource settings

Done when:

- lint is clean or close to clean
- the pipeline structure matches the intended `nf-core` style

## Recommended next concrete implementation step

The best next step after planning is:

- implement Milestone 1

That means editing:

- `nextflow.config`
- `nextflow_schema.json`
- the template input schema file

This is the highest-leverage implementation step because every later module and test depends on the interface being correct.

## Bottom line

The roadmap should stay disciplined:

- define the interface first
- install standard modules second
- wire generic workflow pieces third
- add the local ISAR step after the generic backbone is working
- lock everything down with tests and lint

That sequence gives the project the best chance of becoming a clean and maintainable `nf-core` pipeline instead of a direct copy of the prototype repo.
