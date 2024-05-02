
process COMBINE_SITES {
    cpus 1
    memory '1 GB'
    time '1 h'

    input:
        path(sites)

    output:
        path(out)

    script:
    out = 'sites.combined.gz'
    """
    zcat $sites | sort -t \$'\\t' -k1,1 -k2n,2 | uniq | gzip > $out
    """
}