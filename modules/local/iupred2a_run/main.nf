process IUPRED2A_RUN {
    tag "iupred2a"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container params.iupred2a_container ?: 'docker.io/btrspg/iupred2a:2a'

    input:
    path aa_fasta
    path convert_script

    output:
    path "iupred2a/iupred2a_anchor2.raw", emit: raw
    path "iupred2a/iupred2a_anchor2_isar.out", emit: results
    path "iupred2a/iupred2a.log", emit: log
    tuple val("${task.process}"), val('iupred2a'), val('2a'), emit: versions_iupred2a, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    mkdir -p iupred2a
    fasta="\$(readlink -f ${aa_fasta})"

    if [ ! -s "\$fasta" ]; then
        echo "No amino-acid FASTA records were available for IUPred2A." > iupred2a/iupred2a.log
        touch iupred2a/iupred2a_anchor2.raw
        touch iupred2a/iupred2a_anchor2_isar.out
        exit 0
    fi

    python3 /opt/iupred2a/iupred2a.py -a "\$fasta" long ${args} \\
        > iupred2a/iupred2a_anchor2.raw \\
        2> iupred2a/iupred2a.log

    python3 ${convert_script} \\
        iupred2a/iupred2a_anchor2.raw \\
        "\$fasta" \\
        iupred2a/iupred2a_anchor2_isar.out \\
        >> iupred2a/iupred2a.log 2>&1
    """

    stub:
    """
    mkdir -p iupred2a
    cat > iupred2a/iupred2a_anchor2.raw <<'EOF'
1 M 0.10 0.01
2 A 0.20 0.02
EOF
    cat > iupred2a/iupred2a_anchor2_isar.out <<'EOF'
>stub_isoform
1	M	0.10	0.01
2	A	0.20	0.02
EOF
    echo "stub iupred2a" > iupred2a/iupred2a.log
    """
}
