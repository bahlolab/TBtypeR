
process make_snp_db {
    cpus 2
    memory '4GB'
    time '1h' 
    container 'bahlolab/tbmixbench:quanttb'

    input:
        path(snps)

    output:
        path("output/TBtypeR_QuantTB_snp.db")

    script:
    """
    export PYTHON_EGG_CACHE="\$PWD/python_egg_cache"
    quanttb makesnpdb -g $snps -o 'TBtypeR_QuantTB_snp'
    """
}
