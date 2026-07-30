rule genomad_taxonomy:
    input:
        filtered = out("{sample}","intermediate","final_phage_sequences.fasta")
    output:
        out("{sample}","final_results", "taxonomy","genomad", "final_phage_sequences_summary", "final_phage_sequences_virus_summary.tsv")
    params:
        outdir = out("{sample}","final_results", "taxonomy", "genomad"),
        db = "databases/genomad_db" 
    conda:
         "genomad_env"
    threads: config["resources"]["threads"]
    shell:
        """
        genomad end-to-end --cleanup {input.filtered} {params.outdir} {params.db} --threads {threads}
        """