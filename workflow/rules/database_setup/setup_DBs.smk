rule all:
    input:
        "databases/virsorter2/.complete",
        "databases/dvf_model/.complete",
        "databases/ct3_DBs/.complete",
        "databases/checkv-db/.complete",
        "databases/genomad_db/.complete"

rule download_virsorter2db:
    output:
        touch("databases/virsorter2/.complete")
    conda:
        "envs/virsorter2.yaml"
    shell:
        """
        virsorter setup -d databases/virsorter2 -j 4 
        touch {output}
        """

rule download_dvf_model:
    output:
        touch("databases/dvf_model/.complete")
    shell:
        """
        git clone https://github.com/jessieren/DeepVirFinder.git databases/dvf_model
        touch {output}
        """

rule download_cnt3db:
    output:
        touch("databases/ct3_DBs/.complete")
    conda:
        "envs/cenote-taker3.yaml"
    shell:
        """
        get_ct3_dbs -o databases/ct3_DBs --hmm T --hallmark_tax T --refseq_tax T --mmseqs_cdd T --domain_list T
        touch {output}
        """

rule download_checkvdb:
    output:
        touch("databases/checkv-db/.complete")
    conda:
        "envs/checkv.yaml"
    shell:
        """
        checkv download_database databases/checkv-db
        touch {output}
        """

rule download_genomaddb:
    output:
        touch("databases/genomad_db/.complete")
    conda:
        "envs/genomad.yaml"
    shell:
        """
        genomad download-database databases/genomad_db
        touch {output}
        """
