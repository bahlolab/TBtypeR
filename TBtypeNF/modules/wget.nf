
process WGET {
    cpus 1
    memory '1 GB'
    time '4 h'
    tag "$meta"

    input:
    tuple val(meta), val(url)

    output:
    tuple val(meta), path(output)

    script:
    output = file(url).name
    """
    wget $url
    """
}