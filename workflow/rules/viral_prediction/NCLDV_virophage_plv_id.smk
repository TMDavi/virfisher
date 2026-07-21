rule BEREN:
    input:
        scaffolds=out("{sample}", "intermediate","metaspades", "scaffolds.fasta"),
        metabat_coverage=out("{sample}", "final_results", "read_alignment","{sample}_metabat.coverage")
    
    output:
        summary=out("{sample}", "intermediate","BEREN","Final_results","Run_Summary.txt")

    params:
        out_dir=out("{sample}", "intermediate","BEREN")
    
    conda:
        "BEREN"
    
    threads: config["resources"]["threads"]
    shell:
        """
        cd programs/BEREN

        python BEREN.py -i ../../{input.scaffolds} -o ../../{params.out_dir} -m all -cov ../../{input.metabat_coverage} -t {threads}
        """