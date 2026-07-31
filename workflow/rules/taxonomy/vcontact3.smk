rule vcontact3:
    input:
        out("merged_results","dereplicated","derep_miuvigs.fasta")
    output:
        classification = out("merged_results","vcontact3","exports","final_assignments.csv")
    params:
        outdir = out("merged_results","vcontact3"),
        db = "databases/vcontact3_db" 
    conda:
         "vcontact3"
    threads: config["resources"]["threads"]
    shell:
        """
        vcontact3 run --nucleotide {input} --output {params.outdir} -t {threads} --db-path {params.db}

        rm {params.outdir}/*.faa
        rm {params.outdir}/*.parquet
        rm {params.outdir}/*.gz
        rm {params.outdir}/*.h5

        """