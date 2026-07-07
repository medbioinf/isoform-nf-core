process ISOFORM_GO_ENRICHMENT {
    tag "go_enrichment"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/jungwooseok/webgestalt:1.0.3'

    input:
    path go_script
    path isar_results
    path reference_gene_file

    output:
    path "go_input"         , emit: input
    path "go_enrichment.csv", emit: enrichment
    path "go_summary.txt"   , emit: summary
    path "webgestalt_report", emit: report, optional: true
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def reference_arg = reference_gene_file ? "--reference-gene-file ${reference_gene_file}" : ''
    """
    Rscript ${go_script} \\
        --gene-score-file ${isar_results}/go_gene_scores.csv \\
        --outdir . \\
        --organism ${params.go_organism} \\
        --gene-id-type ${params.go_gene_id_type} \\
        --database ${params.go_database} \\
        --qvalue-cutoff ${params.isar_qvalue_cutoff} \\
        --dif-cutoff ${params.isar_dif_cutoff} \\
        ${reference_arg} \\
        ${args}
    """

    stub:
    """
    mkdir -p go_input
    touch go_input/significant_genes.txt
    touch go_input/background_genes.txt
    cat <<'END_SCORES' > go_input/go_gene_scores.csv
    gene,gene_id,gene_name,min_isoform_switch_q_value,max_abs_dIF,significant_isoform_switch,n_isoforms_tested,source_qvalue_column,source_dif_column
    END_SCORES
    cat <<'END_ENRICHMENT' > go_enrichment.csv
    geneSet,description,size,overlap,enrichmentRatio,pValue,FDR,overlapId,database
    END_ENRICHMENT
    cat <<'END_SUMMARY' > go_summary.txt
    GO enrichment summary
    Status: stub
    END_SUMMARY
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: stub
        webgestaltr: stub
    END_VERSIONS
    """
}
