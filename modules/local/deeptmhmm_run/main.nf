process DEEPTMHMM_RUN {
    tag "deeptmhmm"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/deeptmhmm/deeptmhmm:latest'

    input:
    path aa_fasta
    path convert_script

    output:
    path "deeptmhmm/deeptmhmm_regions_isar.tsv", emit: results
    path "deeptmhmm/predicted_topologies.3line", emit: raw
    path "deeptmhmm/deeptmhmm.log", emit: log
    tuple val("${task.process}"), val('deeptmhmm'), val('latest'), emit: versions_deeptmhmm, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    mkdir -p deeptmhmm
    fasta="\$(readlink -f ${aa_fasta})"

    if [ ! -s "\$fasta" ]; then
        echo "No amino-acid FASTA records were available for DeepTMHMM." > deeptmhmm/deeptmhmm.log
        touch deeptmhmm/predicted_topologies.3line deeptmhmm/deeptmhmm_regions_isar.tsv
        exit 0
    fi

    mkdir -p deeptmhmm_work
    find /openprotein -mindepth 1 -maxdepth 1 ! -name predict.py -exec ln -sfn {} deeptmhmm_work/ \\;
    cp /openprotein/predict.py deeptmhmm_work/predict.py
    perl -pi -e "s#/deeptmhmm_results\\.md#deeptmhmm_results.md#g" deeptmhmm_work/predict.py
    (
        cd deeptmhmm_work
        python3 predict.py --fasta "\$fasta" ${args}
    ) > deeptmhmm/deeptmhmm.log 2>&1
    cp deeptmhmm_work/predicted_topologies.3line deeptmhmm/predicted_topologies.3line

    python3 ${convert_script} \\
        deeptmhmm/predicted_topologies.3line \\
        deeptmhmm/deeptmhmm_regions_isar.tsv
    """

    stub:
    """
    mkdir -p deeptmhmm
    cat > deeptmhmm/predicted_topologies.3line <<'EOF'
>stub_isoform
MKTLLAILAVATLALA
SSSSSSSSSSOOOOOO
EOF
    cat > deeptmhmm/deeptmhmm_regions_isar.tsv <<'EOF'
stub_isoform	signal	1	10
stub_isoform	outside	11	16
EOF
    echo "stub deeptmhmm" > deeptmhmm/deeptmhmm.log
    """
}
