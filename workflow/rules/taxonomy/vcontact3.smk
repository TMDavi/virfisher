rule vcontact3:
    input:
        out("dereplicated","derep_miuvigs.fasta")
    output:
        classification = out("merged_results","vcontact3","final_assignments.csv")
    params:
        outdir = out("vcontact3"),
        finaldir = out("merged_results","vcontact3"),
        db = "databases/vcontact3_db" 
    conda:
         "vcontact3"
    threads: config["resources"]["threads"]
    shell:
        """
        vcontact3 run --nucleotide {input} --output {params.outdir} -t {threads} --db-path {params.db}

        cp {params.outdir}/export/*.csv {params.finaldir}

        rm -r {params.outdir}
        """