configfile: "config.yaml"
import os

rule all:
    input:
        expand("results/{sample}/prodigal/{sample}_proteins.faa", sample=config["samples"].keys())

rule prodigal_gv:
    input:
        assembly = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample])
    output:
        "results/{sample}/prodigal/{sample}_proteins.faa"
    
    conda:
        "/MP_Data/mambaforge/envs/ncdlv_msearch"
    threads: 40    
    shell:
        """
        python scripts/parallel-prodigal-gv.py -t {threads} -q -i {input.assembly} -a {output}
        """