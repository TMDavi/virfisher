rule pharokka:
    input:
        out("{sample}", "intermediate", "final_phage_sequences.fasta")
    output:
        summary = out("{sample}","final_results", "annotation","pharokka","pharokka_cds_final_merged_output.tsv")
    params:
        db = "databases/pharokka_db", 
        outdir = out("{sample}","final_results", "annotation")
    conda:
         "pharokka"
    threads: config["resources"]["threads"]
    shell:
        """
        mkdir -p {params.outdir}
        pharokka run -i {input} -o {params.outdir}/pharokka/ -d {params.db} -t {threads} -f
        """