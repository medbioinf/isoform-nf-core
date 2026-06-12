process PFAM_PREPARE {
    tag "pfam_prepare"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' :
        'quay.io/biocontainers/bioconductor-isoformswitchanalyzer:2.6.0--r44h3df3fcb_0' }"

    input:
    path prepare_script
    path isar_results

    output:
    path "pfam_prepare", emit: results
    path "pfam_prepare/*_AA.fasta", optional: true, emit: aa_fasta
    tuple val("${task.process}"), val('bioconductor-isoformswitchanalyzer'), eval("Rscript -e \"cat(as.character(packageVersion('IsoformSwitchAnalyzeR')))\""), emit: versions_isar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${prepare_script} \\
        --isar-dir ${isar_results} \\
        --outdir pfam_prepare \\
        --top-n ${params.pfam_top_n} \\
        --qvalue-cutoff ${params.isar_qvalue_cutoff} \\
        --dif-cutoff ${params.isar_dif_cutoff} \\
        ${args}
    """

    stub:
    """
    mkdir -p pfam_prepare
    touch pfam_prepare/pfam_prepare_notes.txt
    touch pfam_prepare/pfam_top_candidates.csv
    printf ">stub_isoform\\nMA\\n" > pfam_prepare/isoform_pfam_candidates_AA.fasta
    """
}
