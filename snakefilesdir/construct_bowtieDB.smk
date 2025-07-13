configfile: "config.yaml"
import os

rule all:
    input:
        "bowtieDB/contigs.1.bt2l"

rule rename_contigs:
    input:
        assembly = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample])
    output:
        temp("results/{sample}/{sample}_scaffolds.fasta")
    params:
        sample_name = lambda wildcards: wildcards.sample
    shell:
        "python scripts/rename_contigs.py -i {input.assembly} -o {output} -p {params.sample_name}"

rule concatenate_contigs:
    input:
        expand("results/{sample}/{sample}_scaffolds.fasta", sample=list(config["samples"].keys()))
    output:
        "results/concatenated_scaffolds.fasta"
    shell:
        "cat {input} > {output}"

rule construct_bowtieDB:
    input:
        "results/concatenated_scaffolds.fasta"
    output:
        "bowtieDB/contigs.1.bt2l"
   params:
        db_name = "bowtieDB/contigs"
    threads: 30
    shell:
        "bowtie2-build -f {input} {params.db_name} --large-index --threads {threads}"