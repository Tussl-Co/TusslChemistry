process TCHEM_DETECT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/tussl-chemistry:0.1.0--FILL_IN'
        : 'quay.io/biocontainers/tussl-chemistry:0.1.0--FILL_IN'}"

    input:
    tuple val(meta), path(reads)
    path whitelists

    output:
    tuple val(meta), path("*.chemistry.json"), emit: json
    tuple val("${task.process}"), val('tchem'), eval("tchem --version | sed 's/^tchem //'"), emit: versions_tchem, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // No --emit here: it exits 3 when no command fits, which would fail the task. The report
    // always exits 0 once the reads were read, whatever the verdict.
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tchem \\
        --whitelists ${whitelists} \\
        --format json \\
        ${args} \\
        ${reads} \\
        > ${prefix}.chemistry.json
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    touch ${prefix}.chemistry.json
    """
}
