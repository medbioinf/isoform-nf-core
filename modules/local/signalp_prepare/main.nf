process SIGNALP_PREPARE {
    tag "signalp_prepare"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path prepare_script
    path isar_results

    output:
    path "signalp_prepare", emit: results
    path "signalp_prepare/*_AA.fasta", optional: true, emit: aa_fasta
    tuple val("${task.process}"), val('bioconductor-isoformswitchanalyzer'), eval("Rscript -e \"cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))\""), emit: versions_isar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? " \\\n        ${task.ext.args}" : ''
    """
    Rscript ${prepare_script} \\
        --isar-dir ${isar_results} \\
        --outdir signalp_prepare \\
        --top-n ${params.signalp_top_n} \\
        --qvalue-cutoff ${params.isar_qvalue_cutoff} \\
        --dif-cutoff ${params.isar_dif_cutoff} \\
        --target-genes "${params.annotated_switch_genes ?: ''}"${args}
    """

    stub:
    """
    mkdir -p signalp_prepare
    touch signalp_prepare/signalp_prepare_notes.txt
    touch signalp_prepare/signalp_top_candidates.csv
    printf ">stub_isoform\\nMKTLLA\\n" > signalp_prepare/isoform_signalp_candidates_AA.fasta
    """
}
