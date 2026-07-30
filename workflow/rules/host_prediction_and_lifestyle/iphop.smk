rule iphop:
    input:
        contigs = out("{sample}","intermediate", "final_phage_sequences.fasta")
    params:
        outdir = out("{sample}","final_results","Host_prediction"),
        db_dir = config["iphop_db_dir"] #Adicionar no comando do config
    output:
        out("{sample}","final_results","Host_prediction","Host_prediction_to_genome_m90.csv")
    conda: 
        "iphop_env"
    threads: config["resources"]["threads"]
    shell:
        """
        iphop predict --fa_file {input.contigs} --db_dir {params.db_dir} --out_dir {params.outdir} --num_threads {threads}
        """