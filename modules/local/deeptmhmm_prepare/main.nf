process DEEPTMHMM_PREPARE {
    tag "deeptmhmm_prepare"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path prepare_script
    path isar_results

    output:
    path "deeptmhmm_prepare", emit: results
    path "deeptmhmm_prepare/*_AA.fasta", optional: true, emit: aa_fasta
    tuple val("${task.process}"), val('bioconductor-isoformswitchanalyzer'), eval("Rscript -e \"cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))\""), emit: versions_isar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? " \\\n        ${task.ext.args}" : ''
    """
    Rscript ${prepare_script} \\
        --isar-dir ${isar_results} \\
        --outdir deeptmhmm_prepare \\
        --top-n ${params.deeptmhmm_top_n} \\
        --qvalue-cutoff ${params.isar_qvalue_cutoff} \\
        --dif-cutoff ${params.isar_dif_cutoff} \\
        --target-genes "${params.annotated_switch_genes ?: ''}"${args}
    """

    stub:
    """
    mkdir -p deeptmhmm_prepare
    touch deeptmhmm_prepare/deeptmhmm_prepare_notes.txt
    touch deeptmhmm_prepare/deeptmhmm_top_candidates.csv
    printf ">stub_isoform\\nMKTLLAILAVATLALA\\n" > deeptmhmm_prepare/isoform_deeptmhmm_candidates_AA.fasta
    """
}
