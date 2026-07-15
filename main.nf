#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Anton-Bch/isoform-nf-core
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Anton-Bch/isoform-nf-core
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ISOFORM  } from './workflows/isoform'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_isoform_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_isoform_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow ISOFORM_NF_CORE {

    take:
    samplesheet   // channel: samplesheet read in from --input
    sra_manifest  // channel: SRA manifest read in from --sra_manifest
    metadata_file // path: original metadata CSV used by ISAR

    main:

    //
    // WORKFLOW: Run pipeline
    //
    ISOFORM (
        samplesheet,
        sra_manifest,
        metadata_file
    )
    emit:
    quant_results  = ISOFORM.out.quant_results
    isar_results   = ISOFORM.out.isar_results
    isar_visualization_results = ISOFORM.out.isar_visualization_results
    isar_contrast_summary_results = ISOFORM.out.isar_contrast_summary_results
    go_enrichment_results = ISOFORM.out.go_enrichment_results
    pfam_prepare_results = ISOFORM.out.pfam_prepare_results
    pfam_scan_results = ISOFORM.out.pfam_scan_results
    pfam_import_results = ISOFORM.out.pfam_import_results
    pfam_visualization_results = ISOFORM.out.pfam_visualization_results
    iupred2a_prepare_results = ISOFORM.out.iupred2a_prepare_results
    iupred2a_run_results = ISOFORM.out.iupred2a_run_results
    iupred2a_import_results = ISOFORM.out.iupred2a_import_results
    signalp_prepare_results = ISOFORM.out.signalp_prepare_results
    signalp_run_results = ISOFORM.out.signalp_run_results
    signalp_import_results = ISOFORM.out.signalp_import_results
    deeptmhmm_prepare_results = ISOFORM.out.deeptmhmm_prepare_results
    deeptmhmm_run_results = ISOFORM.out.deeptmhmm_run_results
    deeptmhmm_import_results = ISOFORM.out.deeptmhmm_import_results
    deeploc2_prepare_results = ISOFORM.out.deeploc2_prepare_results
    deeploc2_run_results = ISOFORM.out.deeploc2_run_results
    deeploc2_import_results = ISOFORM.out.deeploc2_import_results
    annotated_switch_plot_results = ISOFORM.out.annotated_switch_plot_results
    multiqc_report = ISOFORM.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions = ISOFORM.out.versions
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.sra_manifest,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    ISOFORM_NF_CORE (
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.sra_manifest,
        PIPELINE_INITIALISATION.out.metadata_file
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        ISOFORM_NF_CORE.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
