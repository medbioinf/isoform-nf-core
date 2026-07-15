//
// Subworkflow with functionality specific to the Anton-Bch/isoform-nf-core pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    sra_manifest      //  string: Path to SRA manifest
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        "",
        "",
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    input_rows = params.input ? samplesheetToList(params.input, "${projectDir}/assets/schema_input.json") : []
    sra_rows = params.sra_manifest ? samplesheetToList(params.sra_manifest, "${projectDir}/assets/schema_sra_manifest.json") : []
    contrast_rows = params.contrasts ? samplesheetToList(params.contrasts, "${projectDir}/assets/schema_contrasts.json") : []
    validateContrasts(input_rows ?: sra_rows, contrast_rows)

    channel
        .fromList(input_rows)
        .map {
            meta, fastq_1, fastq_2 ->
                sampleRowToInput(meta, fastq_1, fastq_2)
        }
        .groupTuple()
        .map { samplesheet ->
            validateInputSamplesheet(samplesheet)
        }
        .map {
            meta, fastqs ->
                return [ meta, fastqs.flatten() ]
        }
        .set { ch_samplesheet }

    channel
        .fromList(sra_rows)
        .map {
            meta, run_accession ->
                sraRowToInput(meta, run_accession)
        }
        .groupTuple()
        .map { manifest ->
            validateInputSamplesheet(manifest)
        }
        .flatMap {
            meta, run_accessions ->
                run_accessions.collect { run_accession ->
                    [ meta + [ sra_run: run_accession ], run_accession ]
                }
        }
        .set { ch_sra_manifest }

    emit:
    samplesheet   = ch_samplesheet
    sra_manifest  = ch_sra_manifest
    metadata_file = channel.value(file(params.input ?: params.sra_manifest, checkIfExists: true))
    versions      = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
    if (params.input && params.sra_manifest) {
        error("Please provide either --input or --sra_manifest, not both")
    }
    if (!params.input && !params.sra_manifest) {
        error("Please provide either --input FASTQ samplesheet or --sra_manifest")
    }
    if (!params.transcript_fasta) {
        error("Please provide --transcript_fasta for Salmon indexing and downstream isoform analysis")
    }
    if (!params.gtf) {
        error("Please provide --gtf so Salmon can produce gene-level mappings and ISAR can map transcripts to genes")
    }

    if (!params.run_isar) {
        def dependent_flags = [
            'run_isar_visualization'    : params.run_isar_visualization,
            'run_isar_contrast_summary' : params.run_isar_contrast_summary,
            'run_go_enrichment'          : params.run_go_enrichment,
            'run_pfam_prepare'           : params.run_pfam_prepare,
            'run_pfam_visualization'     : params.run_pfam_visualization,
            'run_iupred2a'               : params.run_iupred2a,
            'run_signalp'                : params.run_signalp,
            'run_deeptmhmm'              : params.run_deeptmhmm,
            'run_deeploc2'               : params.run_deeploc2,
            'run_annotated_switch_plots' : params.run_annotated_switch_plots,
        ].findAll { name, enabled -> enabled }
        def dependent_inputs = [
            'pfam_results' : params.pfam_results,
            'pfam_db'      : params.pfam_db,
        ].findAll { name, value -> value }

        if (dependent_flags || dependent_inputs) {
            def disable_flags = dependent_flags.keySet().collect { name -> "--${name} false" }
            def remove_inputs = dependent_inputs.keySet().collect { name -> "omit --${name}" }
            def instructions = (disable_flags + remove_inputs).join(', ')
            error("--run_isar false is incompatible with enabled ISAR-dependent options. To run Salmon-only preprocessing, also set or remove: ${instructions}")
        }
    }
}

def sraRowToInput(meta, run_accession) {
    if (!meta.id) {
        error("Please check SRA manifest -> Sample name must be provided")
    }
    if (!run_accession) {
        error("Please check SRA manifest -> run_accession must be provided for sample ${meta.id}")
    }
    return [ meta.id, meta + [ single_end:false ], run_accession ]
}

//
// Validate channels from input samplesheet
//
def sampleRowToInput(meta, fastq_1, fastq_2) {
    if (!meta.id) {
        error("Please check input samplesheet -> Sample name must be provided")
    }

    def resolved_fastq_1 = resolveInputFastq(fastq_1)
    def resolved_fastq_2 = resolveInputFastq(fastq_2)

    if (!resolved_fastq_2) {
        return [ meta.id, meta + [ single_end:true ], [ resolved_fastq_1 ] ]
    } else {
        return [ meta.id, meta + [ single_end:false ], [ resolved_fastq_1, resolved_fastq_2 ] ]
    }
}

def resolveInputFastq(path_value) {
    if (!path_value || (path_value instanceof List && path_value.isEmpty())) {
        return null
    }

    def value = path_value.toString()
    if (!value || value == '[]') {
        return null
    }

    def relative_value = value
    def launch_dir = workflow.launchDir?.toString()
    if (launch_dir && value.startsWith("${launch_dir}/")) {
        relative_value = value.substring(launch_dir.size() + 1)
    }

    def candidates = [value]
    if (!new File(relative_value).isAbsolute()) {
        candidates << "${projectDir}/${relative_value}"
        if (params.input) {
            candidates << "${file(params.input).parent}/${relative_value}"
        }
    }

    def resolved_candidate = candidates.unique()
        .collect { candidate -> file(candidate) }
        .find { candidate_path -> candidate_path.exists() }
    if (resolved_candidate) {
        return resolved_candidate
    }

    error("Please check input samplesheet -> FASTQ file does not exist: ${relative_value}")
}

def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    // Check that multiple runs of the same sample have identical sample-level metadata.
    ['condition', 'replicate', 'strandedness', 'batch'].each { field ->
        def values = metas.collect { meta ->
            meta.containsKey(field) ? meta[field] : null
        }.findAll { value ->
            value != null && !(value instanceof List && value.isEmpty()) && value.toString() != ''
        }.unique()

        if (values.size() > 1) {
            error("Please check input samplesheet -> Multiple runs of a sample must have the same '${field}' value: ${metas[0].id}")
        }
    }

    return [ metas[0], fastqs ]
}

def validateContrasts(sample_rows, contrast_rows) {
    if (!contrast_rows) {
        return
    }

    def sample_conditions = sample_rows.collect { row ->
        getSamplesheetRowMeta(row).condition
    }.findAll { condition ->
        condition != null && condition.toString() != ''
    }.unique()

    contrast_rows.each { row ->
        def contrast = getSamplesheetRowMeta(row)
        ['case', 'control'].each { field ->
            def condition = contrast[field]
            if (!sample_conditions.contains(condition)) {
                error("Please check contrast file -> Contrast '${contrast.id}' ${field} condition '${condition}' is not present in the samplesheet condition column")
            }
        }
    }
}

def getSamplesheetRowMeta(row) {
    if (row instanceof Map) {
        return row
    }
    if (row instanceof List && row[0] instanceof Map) {
        return row[0]
    }
    error("Please check input files -> Could not extract row metadata from samplesheet")
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def tools = [
        "FastQC (Andrews 2010)",
        "fastp (Chen et al. 2018)",
        "Salmon (Patro et al. 2017)",
    ]
    if (params.run_isar) tools << "IsoformSwitchAnalyzeR (Vitting-Seerup and Sandelin 2019)"
    if (params.run_go_enrichment) tools << "WebGestaltR (Liao et al. 2019)"
    if (params.run_pfam_prepare || params.pfam_results || params.pfam_db) tools << "Pfam (Mistry et al. 2021)"
    if (params.run_iupred2a) tools << "IUPred2A (Meszaros et al. 2018)"
    if (params.run_signalp) tools << "SignalP 5.0 (Almagro Armenteros et al. 2019)"
    if (params.run_deeptmhmm) tools << "DeepTMHMM (Hallgren et al. 2022)"
    if (params.run_deeploc2) tools << "DeepLoc 2.0 (Thumuluri et al. 2022)"
    tools << "MultiQC (Ewels et al. 2016)"
    return "Tools used in the workflow included: ${tools.join(', ')}."
}

def toolBibliographyText() {
    def references = [
        "<li>Andrews S. (2010) FastQC. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.</li>",
        "<li>Chen S, Zhou Y, Chen Y, Gu J. (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. doi: 10.1093/bioinformatics/bty560.</li>",
        "<li>Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. (2017) Salmon provides fast and bias-aware quantification of transcript expression. Nat Methods. doi: 10.1038/nmeth.4197.</li>",
    ]
    if (params.run_isar) references << "<li>Vitting-Seerup K, Sandelin A. (2019) IsoformSwitchAnalyzeR: analysis of changes in genome-wide patterns of alternative isoform usage and its functional consequences. Mol Cancer Res. doi: 10.1158/1541-7786.MCR-18-0262.</li>"
    if (params.run_go_enrichment) references << "<li>Liao Y, Wang J, Jaehnig EJ, Shi Z, Zhang B. (2019) WebGestalt 2019: gene set analysis toolkit with revamped UIs and APIs. Nucleic Acids Res. doi: 10.1093/nar/gkz401.</li>"
    if (params.run_pfam_prepare || params.pfam_results || params.pfam_db) references << "<li>Mistry J, Chuguransky S, Williams L, et al. (2021) Pfam: The protein families database in 2021. Nucleic Acids Res. doi: 10.1093/nar/gkaa913.</li>"
    if (params.run_iupred2a) references << "<li>Meszaros B, Erdos G, Dosztanyi Z. (2018) IUPred2A: context-dependent prediction of protein disorder. Nucleic Acids Res. doi: 10.1093/nar/gky384.</li>"
    if (params.run_signalp) references << "<li>Almagro Armenteros JJ, et al. (2019) SignalP 5.0 improves signal peptide predictions using deep neural networks. Nat Biotechnol. doi: 10.1038/s41587-019-0036-z.</li>"
    if (params.run_deeptmhmm) references << "<li>Hallgren J, et al. (2022) DeepTMHMM predicts alpha and beta transmembrane proteins using deep neural networks. bioRxiv. doi: 10.1101/2022.04.08.487609.</li>"
    if (params.run_deeploc2) references << "<li>Thumuluri V, et al. (2022) DeepLoc 2.0: multi-label subcellular localization prediction using protein language models. Nucleic Acids Res. doi: 10.1093/nar/gkac278.</li>"
    references << "<li>Ewels P, Magnusson M, Lundin S, Kaller M. (2016) MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. doi: 10.1093/bioinformatics/btw354.</li>"
    return references.join(' ').trim()
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
