WORKDIR = config["WORKDIR"]

rule rename_contigs:
    input:  
        assembly=out("{sample}", "intermediate", "metaspades", "scaffolds.fasta")
    output:
        temp(out("{sample}", "intermediate", "metaspades", "{sample}_scaffolds.fasta"))
    params:
        sample_name=lambda wc: wc.sample
    shell:
        """
        python {WORKDIR}/scripts/rename_contigs.py -i {input.assembly} -o {output} -p {params.sample_name}
        """

rule concatenate_contigs:
    input:
        expand(
            out("{sample}", "intermediate", "metaspades", "{sample}_scaffolds.fasta"),
            sample=list(config["samples"].keys())
        )
    output:
        temp(out("mappingDB", "concatenated_scaffolds.fasta"))
    shell:
        """
        cat {input} > {output}
        """

rule filter_contigs_length:
    input:
        out("mappingDB", "concatenated_scaffolds.fasta")
    output:
        out("mappingDB", "scaffolds_1000bp.fasta")
    shell:
        """
        seqkit seq -m 1000 {input} > {output}
        """

rule construct_bwaDB:
    input:
        fasta=out("mappingDB", "scaffolds_1000bp.fasta")
    output:
        expand(out("mappingDB", "mappingDB.{suffix}"),suffix=["amb","ann","bwt","pac","sa"])
    params:
        index_prefix=out("mappingDB", "mappingDB")
    shell:
        """
        bwa index -p {params.index_prefix} {input.fasta} -a bwtsw
        """