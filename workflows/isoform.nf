/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { SALMON_INDEX           } from '../modules/nf-core/salmon/index/main'
include { SALMON_QUANT           } from '../modules/nf-core/salmon/quant/main'
include { FETCH_SRA_FASTQ        } from '../modules/local/fetch_sra_fastq/main'
include { ISAR_ANALYSIS          } from '../modules/local/isar_analysis/main'
include { ISAR_VISUALIZATION     } from '../modules/local/isar_visualization/main'
include { ISAR_CONTRAST_SUMMARY  } from '../modules/local/isar_contrast_summary/main'
include { ISOFORM_GO_ENRICHMENT  } from '../modules/local/isoform_go_enrichment/main'
include { PFAM_PREPARE           } from '../modules/local/pfam_prepare/main'
include { PFAM_SCAN              } from '../modules/local/pfam_scan/main'
include { PFAM_IMPORT            } from '../modules/local/pfam_import/main'
include { PFAM_VISUALIZATION     } from '../modules/local/pfam_visualization/main'
include { IUPRED2A_PREPARE       } from '../modules/local/iupred2a_prepare/main'
include { IUPRED2A_RUN           } from '../modules/local/iupred2a_run/main'
include { IUPRED2A_IMPORT        } from '../modules/local/iupred2a_import/main'
include { SIGNALP_PREPARE        } from '../modules/local/signalp_prepare/main'
include { SIGNALP_RUN            } from '../modules/local/signalp_run/main'
include { SIGNALP_IMPORT         } from '../modules/local/signalp_import/main'
include { DEEPTMHMM_PREPARE      } from '../modules/local/deeptmhmm_prepare/main'
include { DEEPTMHMM_RUN          } from '../modules/local/deeptmhmm_run/main'
include { DEEPTMHMM_IMPORT       } from '../modules/local/deeptmhmm_import/main'
include { DEEPLOC2_PREPARE       } from '../modules/local/deeploc2_prepare/main'
include { DEEPLOC2_RUN           } from '../modules/local/deeploc2_run/main'
include { DEEPLOC2_IMPORT        } from '../modules/local/deeploc2_import/main'
include { ANNOTATED_SWITCH_PLOTS } from '../modules/local/annotated_switch_plots/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_isoform_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ISOFORM {

    take:
    ch_fastq_samplesheet // channel: samplesheet read in from --input
    ch_sra_manifest      // channel: SRA manifest read in from --sra_manifest
    ch_metadata_file     // path: original metadata CSV used by ISAR
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // MODULE: Optionally download SRA runs and normalize them to FASTQ tuples
    //
    FETCH_SRA_FASTQ (
        ch_sra_manifest
    )
    ch_sra_reads = FETCH_SRA_FASTQ.out.reads
        .map { meta, reads ->
            def ordered_reads = (reads instanceof List ? reads : [reads]).sort { it.name }
            if (!(ordered_reads.size() in [1, 2])) {
                error("SRA run ${meta.sra_run} for sample ${meta.id} produced ${ordered_reads.size()} FASTQ files; expected one single-end file or one paired-end file pair")
            }
            [ meta.id, [ meta, ordered_reads ] ]
        }
        .groupTuple()
        .map { sample_id, run_records ->
            def sorted_runs = run_records.sort { left, right ->
                left[0].sra_run.toString() <=> right[0].sra_run.toString()
            }
            def read_counts = sorted_runs.collect { run_record -> run_record[1].size() }.unique()
            if (read_counts.size() != 1) {
                error("SRA runs for sample ${sample_id} mix single-end and paired-end data; all runs for a biological sample must have the same endedness")
            }

            def sample_meta = new LinkedHashMap(sorted_runs[0][0])
            sample_meta.remove('sra_run')
            sample_meta.single_end = read_counts[0] == 1

            // CAT_FASTQ separates paired reads by alternating list position, so retain
            // [run1_R1, run1_R2, run2_R1, run2_R2, ...] ordering after sorting runs.
            def sample_reads = sorted_runs.collectMany { run_record -> run_record[1] }
            [ sample_meta, sample_reads ]
        }
    ch_samplesheet = ch_fastq_samplesheet.mix(ch_sra_reads)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Concatenate multiple sequencing runs per sample
    //
    CAT_FASTQ (
        ch_samplesheet
    )

    //
    // MODULE: Trim adapters and low-quality sequence
    //
    FASTP (
        CAT_FASTQ.out.reads.map { meta, reads -> [ meta, reads, [] ] },
        false,
        false,
        false
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.log.collect{it[1]})

    //
    // MODULE: Build Salmon transcriptome index
    //
    ch_genome_fasta = params.genome_fasta ?
        channel.value(file(params.genome_fasta, checkIfExists: true)) :
        channel.value([])
    ch_transcript_fasta = channel.value(file(params.transcript_fasta, checkIfExists: true))

    SALMON_INDEX (
        ch_genome_fasta,
        ch_transcript_fasta
    )

    //
    // MODULE: Quantify transcript abundance per sample
    //
    ch_gtf = channel.value(file(params.gtf, checkIfExists: true))
    SALMON_QUANT (
        FASTP.out.reads,
        SALMON_INDEX.out.index,
        ch_gtf,
        ch_transcript_fasta,
        false,
        params.salmon_lib_type ?: ''
    )
    ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.json_info.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.lib_format_counts.collect{it[1]})

    //
    // MODULE: Import Salmon quantifications into IsoformSwitchAnalyzeR
    //
    ch_isar_results = channel.empty()
    ch_isar_visualization_results = channel.empty()
    ch_isar_contrast_summary_results = channel.empty()
    ch_go_enrichment_results = channel.empty()
    ch_pfam_prepare_results = channel.empty()
    ch_pfam_scan_results = channel.empty()
    ch_pfam_import_results = channel.empty()
    ch_pfam_visualization_results = channel.empty()
    ch_iupred2a_prepare_results = channel.empty()
    ch_iupred2a_run_results = channel.empty()
    ch_iupred2a_import_results = channel.empty()
    ch_signalp_prepare_results = channel.empty()
    ch_signalp_run_results = channel.empty()
    ch_signalp_import_results = channel.empty()
    ch_deeptmhmm_prepare_results = channel.empty()
    ch_deeptmhmm_run_results = channel.empty()
    ch_deeptmhmm_import_results = channel.empty()
    ch_deeploc2_prepare_results = channel.empty()
    ch_deeploc2_run_results = channel.empty()
    ch_deeploc2_import_results = channel.empty()
    ch_annotated_switch_plot_results = channel.empty()
    if (params.run_isar) {
        ch_analysis_script = channel.value(file("${projectDir}/bin/run_isar_analysis.R", checkIfExists: true))
        ch_input_samplesheet = ch_metadata_file
        ch_contrasts = params.contrasts ?
            channel.value(file(params.contrasts, checkIfExists: true)) :
            channel.value([])
        ch_quant_dirs = SALMON_QUANT.out.results.map { meta, quant_dir -> quant_dir }.collect()

        ISAR_ANALYSIS (
            ch_analysis_script,
            ch_input_samplesheet,
            ch_contrasts,
            ch_quant_dirs,
            ch_gtf,
            ch_transcript_fasta
        )
        ch_isar_results = ISAR_ANALYSIS.out.results
        ch_current_annotated_isar = ISAR_ANALYSIS.out.results

        if (params.run_isar_visualization) {
            ch_visualization_script = channel.value(file("${projectDir}/bin/run_isar_visualization.R", checkIfExists: true))
            ISAR_VISUALIZATION (
                ch_visualization_script,
                ISAR_ANALYSIS.out.results
            )
            ch_isar_visualization_results = ISAR_VISUALIZATION.out.results
        }

        if (params.run_isar_contrast_summary) {
            ch_contrast_summary_script = channel.value(file("${projectDir}/bin/run_isar_contrast_summary.R", checkIfExists: true))
            ISAR_CONTRAST_SUMMARY (
                ch_contrast_summary_script,
                ISAR_ANALYSIS.out.results
            )
            ch_isar_contrast_summary_results = ISAR_CONTRAST_SUMMARY.out.results
        }

        if (params.run_go_enrichment) {
            ch_go_enrichment_script = channel.value(file("${projectDir}/bin/run_go_enrichment.R", checkIfExists: true))
            ch_go_reference_gene_file = params.go_reference_gene_file ?
                channel.value(file(params.go_reference_gene_file, checkIfExists: true)) :
                channel.value([])
            ISOFORM_GO_ENRICHMENT (
                ch_go_enrichment_script,
                ISAR_ANALYSIS.out.results,
                ch_go_reference_gene_file
            )
            ch_go_enrichment_results = ISOFORM_GO_ENRICHMENT.out.enrichment
            ch_multiqc_files = ch_multiqc_files.mix(ISOFORM_GO_ENRICHMENT.out.enrichment)
            ch_multiqc_files = ch_multiqc_files.mix(ISOFORM_GO_ENRICHMENT.out.summary)
            ch_versions = ch_versions.mix(ISOFORM_GO_ENRICHMENT.out.versions)
        }

        if (params.run_pfam_prepare || (!params.pfam_results && params.pfam_db)) {
            ch_pfam_prepare_script = channel.value(file("${projectDir}/bin/run_pfam_prepare.R", checkIfExists: true))
            PFAM_PREPARE (
                ch_pfam_prepare_script,
                ISAR_ANALYSIS.out.results
            )
            ch_pfam_prepare_results = PFAM_PREPARE.out.results
        }

        if (!params.pfam_results && params.pfam_db) {
            ch_pfam_db = channel.value(file(params.pfam_db, checkIfExists: true))
            PFAM_SCAN (
                PFAM_PREPARE.out.aa_fasta,
                ch_pfam_db
            )
            ch_pfam_scan_results = PFAM_SCAN.out.results
        }

        if (params.pfam_results || params.pfam_db) {
            ch_pfam_import_script = channel.value(file("${projectDir}/bin/run_pfam_import.R", checkIfExists: true))
            ch_pfam_results = params.pfam_results ?
                channel.value(file(params.pfam_results, checkIfExists: true)) :
                PFAM_SCAN.out.results
            PFAM_IMPORT (
                ch_pfam_import_script,
                ch_current_annotated_isar,
                ch_pfam_results
            )
            ch_pfam_import_results = PFAM_IMPORT.out.results
            ch_current_annotated_isar = PFAM_IMPORT.out.results

            if (params.run_pfam_visualization) {
                ch_pfam_visualization_script = channel.value(file("${projectDir}/bin/run_pfam_visualization.R", checkIfExists: true))
                PFAM_VISUALIZATION (
                    ch_pfam_visualization_script,
                    PFAM_IMPORT.out.results
                )
                ch_pfam_visualization_results = PFAM_VISUALIZATION.out.results
            }
        }

        if (params.run_iupred2a) {
            ch_iupred2a_prepare_script = channel.value(file("${projectDir}/bin/run_iupred2a_prepare.R", checkIfExists: true))
            ch_iupred2a_convert_script = channel.value(file("${projectDir}/bin/convert_iupred2a_to_isar.py", checkIfExists: true))
            ch_iupred2a_import_script = channel.value(file("${projectDir}/bin/run_iupred2a_import.R", checkIfExists: true))

            IUPRED2A_PREPARE (
                ch_iupred2a_prepare_script,
                ISAR_ANALYSIS.out.results
            )
            ch_iupred2a_prepare_results = IUPRED2A_PREPARE.out.results

            IUPRED2A_RUN (
                IUPRED2A_PREPARE.out.aa_fasta,
                ch_iupred2a_convert_script
            )
            ch_iupred2a_run_results = IUPRED2A_RUN.out.results

            IUPRED2A_IMPORT (
                ch_iupred2a_import_script,
                ch_current_annotated_isar,
                IUPRED2A_RUN.out.results
            )
            ch_iupred2a_import_results = IUPRED2A_IMPORT.out.results
            ch_current_annotated_isar = IUPRED2A_IMPORT.out.results
        }

        if (params.run_signalp) {
            ch_signalp_prepare_script = channel.value(file("${projectDir}/bin/run_signalp_prepare.R", checkIfExists: true))
            ch_signalp_import_script = channel.value(file("${projectDir}/bin/run_signalp_import.R", checkIfExists: true))

            SIGNALP_PREPARE (
                ch_signalp_prepare_script,
                ISAR_ANALYSIS.out.results
            )
            ch_signalp_prepare_results = SIGNALP_PREPARE.out.results

            SIGNALP_RUN (
                SIGNALP_PREPARE.out.aa_fasta
            )
            ch_signalp_run_results = SIGNALP_RUN.out.results

            SIGNALP_IMPORT (
                ch_signalp_import_script,
                ch_current_annotated_isar,
                SIGNALP_RUN.out.results
            )
            ch_signalp_import_results = SIGNALP_IMPORT.out.results
            ch_current_annotated_isar = SIGNALP_IMPORT.out.results
        }

        if (params.run_deeptmhmm) {
            ch_deeptmhmm_prepare_script = channel.value(file("${projectDir}/bin/run_deeptmhmm_prepare.R", checkIfExists: true))
            ch_deeptmhmm_convert_script = channel.value(file("${projectDir}/bin/convert_deeptmhmm_3line_to_isar.py", checkIfExists: true))
            ch_deeptmhmm_import_script = channel.value(file("${projectDir}/bin/run_deeptmhmm_import.R", checkIfExists: true))

            DEEPTMHMM_PREPARE (
                ch_deeptmhmm_prepare_script,
                ISAR_ANALYSIS.out.results
            )
            ch_deeptmhmm_prepare_results = DEEPTMHMM_PREPARE.out.results

            DEEPTMHMM_RUN (
                DEEPTMHMM_PREPARE.out.aa_fasta,
                ch_deeptmhmm_convert_script
            )
            ch_deeptmhmm_run_results = DEEPTMHMM_RUN.out.results

            DEEPTMHMM_IMPORT (
                ch_deeptmhmm_import_script,
                ch_current_annotated_isar,
                DEEPTMHMM_RUN.out.results
            )
            ch_deeptmhmm_import_results = DEEPTMHMM_IMPORT.out.results
            ch_current_annotated_isar = DEEPTMHMM_IMPORT.out.results
        }

        if (params.run_deeploc2) {
            ch_deeploc2_prepare_script = channel.value(file("${projectDir}/bin/run_deeploc2_prepare.R", checkIfExists: true))
            ch_deeploc2_convert_script = channel.value(file("${projectDir}/bin/convert_deeploc2_to_isar.py", checkIfExists: true))
            ch_deeploc2_import_script = channel.value(file("${projectDir}/bin/run_deeploc2_import.R", checkIfExists: true))

            DEEPLOC2_PREPARE (
                ch_deeploc2_prepare_script,
                ISAR_ANALYSIS.out.results
            )
            ch_deeploc2_prepare_results = DEEPLOC2_PREPARE.out.results

            DEEPLOC2_RUN (
                DEEPLOC2_PREPARE.out.aa_fasta,
                ch_deeploc2_convert_script
            )
            ch_deeploc2_run_results = DEEPLOC2_RUN.out.results

            DEEPLOC2_IMPORT (
                ch_deeploc2_import_script,
                ch_current_annotated_isar,
                DEEPLOC2_RUN.out.results
            )
            ch_deeploc2_import_results = DEEPLOC2_IMPORT.out.results
            ch_current_annotated_isar = DEEPLOC2_IMPORT.out.results
        }

        if (params.run_annotated_switch_plots) {
            ch_annotated_switch_plot_script = channel.value(file("${projectDir}/bin/run_annotated_switch_plots.R", checkIfExists: true))
            ANNOTATED_SWITCH_PLOTS (
                ch_annotated_switch_plot_script,
                ch_current_annotated_isar
            )
            ch_annotated_switch_plot_results = ANNOTATED_SWITCH_PLOTS.out.results
        }
    }

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'isoform_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    quant_results  = SALMON_QUANT.out.results       // channel: [ meta, path(salmon_quant_dir) ]
    isar_results   = ch_isar_results                // channel: path(isar_analysis)
    isar_visualization_results = ch_isar_visualization_results // channel: path(isar_visualization)
    isar_contrast_summary_results = ch_isar_contrast_summary_results // channel: path(isar_contrast_summary)
    go_enrichment_results = ch_go_enrichment_results // channel: path(go_enrichment.csv)
    pfam_prepare_results = ch_pfam_prepare_results  // channel: path(pfam_prepare)
    pfam_scan_results = ch_pfam_scan_results        // channel: path(pfam_scan.out)
    pfam_import_results = ch_pfam_import_results    // channel: path(pfam_import)
    pfam_visualization_results = ch_pfam_visualization_results // channel: path(pfam_visualization)
    iupred2a_prepare_results = ch_iupred2a_prepare_results // channel: path(iupred2a_prepare)
    iupred2a_run_results = ch_iupred2a_run_results  // channel: path(iupred2a_anchor2_isar.out)
    iupred2a_import_results = ch_iupred2a_import_results // channel: path(iupred2a_import)
    signalp_prepare_results = ch_signalp_prepare_results // channel: path(signalp_prepare)
    signalp_run_results = ch_signalp_run_results    // channel: path(signalp5_summary.signalp5)
    signalp_import_results = ch_signalp_import_results // channel: path(signalp_import)
    deeptmhmm_prepare_results = ch_deeptmhmm_prepare_results // channel: path(deeptmhmm_prepare)
    deeptmhmm_run_results = ch_deeptmhmm_run_results // channel: path(deeptmhmm_regions_isar.tsv)
    deeptmhmm_import_results = ch_deeptmhmm_import_results // channel: path(deeptmhmm_import)
    deeploc2_prepare_results = ch_deeploc2_prepare_results // channel: path(deeploc2_prepare)
    deeploc2_run_results = ch_deeploc2_run_results // channel: path(deeploc2_isar.csv)
    deeploc2_import_results = ch_deeploc2_import_results // channel: path(deeploc2_import)
    annotated_switch_plot_results = ch_annotated_switch_plot_results // channel: path(annotated_switch_plots)
    multiqc_report = MULTIQC.out.report.toList()    // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
