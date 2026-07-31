rule miuvig_quality:
    input:
        out("per_sample_results", "{sample}", "intermediate", "final_phage_sequences.fasta")
    output:
        summary = out("per_sample_results", "{sample}","final_results", "quality_summary", "quality_summary.tsv")
    params:
        db = "databases/checkv-db-v1.5/", 
        outdir = out("per_sample_results", "{sample}","final_results", "quality_summary")
    conda:
         "viral-id-sop"
    threads: config["resources"]["threads"]
    shell:
        """
        checkv end_to_end {input} {params.outdir} -t {threads} -d {params.db}
        rm -r {params.outdir}/tmp
        """