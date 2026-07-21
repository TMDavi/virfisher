rule prodigal_gv:
    input:
        assembly = out("{sample}", "intermediate","metaspades", "scaffolds.fasta")
    output:
        out("{sample}", "intermediate","prodigal","{sample}_proteins.faa")
    conda:
        "ncdlv_msearch"
    threads: 40    
    shell:
        """
        python {WORKDIR}/scripts/parallel-prodigal-gv.py -t {threads} -q -i {input.assembly} -a {output}
        """