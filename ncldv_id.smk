configfile: "config.yaml"
import os

rule all:
    input:
        expand("results/{sample}/ncldv_id/{sample}.full_output.txt", sample=config["samples"].keys())

rule prodigal_gv:
    input:
        assembly = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample])
    output:
        "results/{sample}/prodigal/{sample}.faa"
    threads: 40
    shell:
        """
        python scripts/parallel-prodigal-gv.py -t {threads} -q -i {input.assembly} -a {output}
        """
rule ncldv_msearch:
    input:
        genes = "results/{sample}/prodigal/{sample}.faa"
    output:
        "results/{sample}/ncldv_id/{sample}.full_output.txt"
    params:
        samplename = 
        outdir = "results/{sample}/ncldv_id/"
    threads: 40
    shell:
        "python viral_HMMs/giantDB/ncldv_markersearch/ncldv_markersearch.py -i {input.genes} -n {params.samplename} -t {threads}"


    