workdir = config['workdir']

rule rename_contigs:
    input:
        assembly = out("{sample}", "intermediate","metaspades", "scaffolds.fasta")
    output:
        temp(out("{sample}", "intermediate","metaspades", "{sample}_scaffolds.fasta"))
    params:
        sample_name = lambda wildcards: wildcards.sample
    shell:
        "python {{workdir}}/scripts/rename_contigs.py -i {input.assembly} -o {output} -p {params.sample_name}"

rule concatenate_contigs:
    input:
        expand(out("{sample}", "intermediate","metaspades", "{sample}_scaffolds.fasta"), sample=list(config["samples"].keys()))
    output:
        temp(out("mappingDB", "concatenated_scaffolds.fasta"))
    shell:
        "cat {input} > {output}"

rule filter_contigs_length:
    input:
        out("mappingDB", "concatenated_scaffolds.fasta")
    output:
        out("mappingDB", "scaffolds_1000bp.fasta")
    shell:
        """
        seqkit seq -m 1000 {input} > {output}
        """

rule construct_bowtieDB:
    input:
        out("mappingDB", "scaffolds_1000bp.fasta")
    output:
        expand(
            "mappingDB/mappingDB.{suffix}",
            suffix=[
                "1.bt2l",
                "2.bt2l",
                "3.bt2l",
                "4.bt2l",
                "rev.1.bt2l",
                "rev.2.bt2l",
            ]
        )
    params:
        db_name = "mappingDB/mappingDB"
    threads: 
        config['resourses']['threads']
    shell:
        "bowtie2-build -f {input} {params.db_name} --large-index --threads {threads}"

