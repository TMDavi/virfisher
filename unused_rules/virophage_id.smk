rule run_virophage_affi:
    input:
        assembly = out("{sample}", "intermediate","metaspades", "scaffolds.fasta")
    output:
        out("{sample}","intermediate","virophage_plv_id", "{sample}_final_affiliation.tsv")
    params:
        outdir = out("{sample}","intermediate","virophage_plv_id")
    conda:
        "virophage_classification"
    shell:
        """
        ./programs/ICTV_VirophageSG/run_virophage_affi.pl -i {input.assembly} -o {params.outdir}
        """

rule get_virophage_and_plv_lists:
    input:
        out("{sample}","intermediate","virophage_plv_id", "{sample}_final_affiliation.tsv")
    output:
        out("{sample}","intermediate","virophage_plv_id","putative_virophage_list.tsv"),
        out("{sample}","intermediate","virophage_plv_id","putative_plv_list.tsv")
    params:
        outdir = out("{sample}","intermediate","virophage_plv_id")
    conda:
        "virophage_classification"
    shell:
        """
        python {WORKDIR}/scripts/virophage_extract.py -i {input} -o {params.outdir} 
        """

rule filter_contigs:
    input:
        ids_virophage = out("{sample}","intermediate","virophage_plv_id","putative_virophage_list.tsv"), 
        ids_plv = out("{sample}","intermediate","virophage_plv_id","putative_plv_list.tsv"), 
        fasta = out("{sample}", "intermediate","metaspades", "scaffolds.fasta")
    output:
        merged_ids = out("{sample}","intermediate","virophage_plv_id","merged_ids.tsv"),
        virophages = out("{sample}","intermediate","virophage_plv_id","putative_virophage_plv_contigs.fasta")

    shell:
        """
        cat {input.ids_virophage} {input.ids_plv} |sort | uniq >> {output.merged_ids}
        seqkit grep -f {output.merged_ids} {input.fasta} >> {output.virophages}
        """

rule checkv_viro:
    input:
        assembly = out("{sample}","intermediate","virophage_plv_id","putative_virophage_plv_contigs.fasta")
    output:
        out("{sample}","intermediate","virophage_plv_id","checkv", "viruses.fna")
    params:
        db = "databases/checkv-db-v1.5/",
        outdir = out("{sample}","intermediate","virophage_plv_id","checkv")
    conda:
         "viral-id-sop"
    threads: config["resources"]["threads"]
    shell:
        """
        checkv end_to_end {input} {params.outdir} -t {threads} -d {params.db}
        """

rule genomad_viro:
    input:
        filtered = out("{sample}","intermediate","virophage_plv_id","checkv", "viruses.fna")
    output:
        out("{sample}","intermediate","virophage_plv_id","genomad","viruses_annotate","viruses_genes.tsv")

    params:
        outdir = out("{sample}","intermediate","virophage_plv_id","genomad"),
        db = "databases/genomad_db" #adcionar no config 
    conda:
         "genomad_env"
    threads: config["resources"]["threads"]
    shell:
        """
        genomad end-to-end --cleanup {input.filtered} {params.outdir} {params.db} --threads {threads}
        """