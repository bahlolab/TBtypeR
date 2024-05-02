
process make_bc_snps {
    cpus 1
    memory '2GB'
    time '1h' 
    container 'bahlolab/tbmixbench:quanttb'

    input:
        path(fasta)

    output:
        path("output/temp/*.snps")

    script:
    """
    export PYTHON_EGG_CACHE="\$PWD/python_egg_cache"
    quanttb makesnpdb -k -g $fasta
    """
}
