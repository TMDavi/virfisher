WORKDIR = "../../../envs"

rule all:
    input:
        "databases/virsorter2/.complete",
        "programs/DeepVirFinder/.complete",
        "databases/ct3_DBs/.complete",
        "databases/checkv-db-v1.5/.complete",
        "databases/genomad_db/.complete",
        "databases/vibrant_db/.complete",
        "databases/pharokka_db/.complete",
        "databases/vcontact3/.complete",
        "databases/iphop_db/.complete",
        "databases/dramv_db/.complete"

rule download_virsorter2db:
    output:
        touch("databases/virsorter2/.complete")
    conda:
        f"{WORKDIR}/viral_id_sop.yaml"
    shell:
        """
        virsorter setup -d databases/virsorter2 -j 4 
        touch {output}
        """

rule download_dvf:
    output:
        touch("programs/DeepVirFinder/.complete")
    conda:
        f"{WORKDIR}/deepvirfinder.yaml"
    shell:
        """
        git clone https://github.com/jessieren/DeepVirFinder.git programs/
        touch {output}
        """

rule download_cnt3db:
    output:
        touch("databases/ct3_DBs/.complete")
    conda:
        f"{WORKDIR}/cenote-taker3.yaml"
    shell:
        """
        get_ct3_dbs -o databases/ct3_DBs --hmm T --hallmark_tax T --refseq_tax T --mmseqs_cdd T --domain_list T
        touch {output}
        """

rule download_checkvdb:
    output:
        touch("databases/checkv-db-v1.5/.complete")
    conda:
        f"{WORKDIR}/viral_id_sop.yaml"
    shell:
        """
        checkv download_database databases/
        touch {output}
        """

rule download_genomaddb:
    output:
        touch("databases/genomad_db/.complete")
    conda:
        f"{WORKDIR}/genomad_env.yaml"
    shell:
        """
        genomad download-database databases/genomad_db
        touch {output}
        """
rule vibrant_db:
    output:
        touch("databases/vibrant_db/.complete")
    conda:
        f"{WORKDIR}/vibrant.yaml"
    shell:
        """
        python3 VIBRANT_setup.py
        touch {output}
        """

rule download_pharokka_db:
    output:
        touch("databases/pharokka_db/.complete")
    conda:
        f"{WORKDIR}/pharokka.yaml"
    shell:
        """
        pharokka install -o databases/pharokka_db
        touch {output}
        """

rule download_vcontact3_db:
    output:
        touch("databases/vcontact3/.complete")
    conda:
        f"{WORKDIR}/vcontact3.yaml"
    shell:
        """
        vcontact3 prepare_databases --get-version latest --set-location databases/vcontact3
        touch {output}
        """

rule download_iphop_db:
    output:
        touch("databases/iphop_db/.complete")
    conda:
        f"{WORKDIR}/iphop.yaml"
    shell:
        """
        mkdir databases/iphop_db/
        iphop download --db_dir databases/iphop_db/ -dbv iPHoP_db_rw_for-test
        touch {output}
        """


rule download_dramv_db:
    output:
        touch("databases/dramv_db/.complete")
    conda:
        f"{WORKDIR}/viral_id_sop.yaml"
    shell:
        """
        DRAM-setup.py prepare_databases --skip_uniref --output_dir databases/dramv_db
        touch {output}
        """
