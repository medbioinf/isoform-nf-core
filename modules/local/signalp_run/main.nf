process SIGNALP_RUN {
    tag "signalp"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/btrspg/signalp:5.0b'

    input:
    path aa_fasta

    output:
    path "signalp/signalp5_summary.signalp5", emit: results
    path "signalp/signalp.log", emit: log
    tuple val("${task.process}"), val('signalp'), val('5.0b'), emit: versions_signalp, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    set -euo pipefail

    mkdir -p signalp
    fasta="\$(readlink -f ${aa_fasta})"

    if [ ! -s "\$fasta" ]; then
        echo "No amino-acid FASTA records were available for SignalP." > signalp/signalp.log
        touch signalp/signalp5_summary.signalp5
        exit 0
    fi

    export PATH="/opt/signalp/bin:\$PATH"
    signalp \\
        -fasta "\$fasta" \\
        -org euk \\
        -format short \\
        -stdout \\
        -prefix signalp/signalp5 \\
        ${args} \\
        > signalp/signalp5_summary.signalp5 \\
        2> signalp/signalp.log
    """

    stub:
    """
    mkdir -p signalp
    cat > signalp/signalp5_summary.signalp5 <<'EOF'
# SignalP-5.0	Organism: Eukarya	Timestamp: stub
# ID	Prediction	SP(Sec/SPI)	OTHER	CS Position
stub_isoform	OTHER	0.0000	1.0000
EOF
    echo "stub signalp" > signalp/signalp.log
    """
}
