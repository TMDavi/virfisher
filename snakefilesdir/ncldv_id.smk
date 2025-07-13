configfile: "config.yaml"
import os

rule all:
    input:
        expand("results/{sample}/ncldv_id/{sample}.full_output.txt", sample=config["samples"].keys())

rule ncldv_msearch:
    input:
        genes = "results/{sample}/prodigal/{sample}_proteins.faa"
    output:
        "results/{sample}/ncldv_id/{sample}.full_output.txt"
    params:
        input_dir = "results/{sample}/prodigal",
        samplename = lambda wildcards: wildcards.sample,
        outdir = "results/{sample}/ncldv_id/"
        
    conda:
        "/MP_Data/mambaforge/envs/ncdlv_msearch"
    threads: 40    
    shell:
        "python viral_HMMs/giantDB/ncldv_markersearch/ncldv_markersearch.py -i {params.input_dir} -n {params.samplename} -t {threads} --allhits -o {params.outdir}"


#Interpretar output
#Filtrr os contigs que são provaveis virus gigantes
#ROdar ViralRecall pra remover contaminações
#Lista final de prováveis vírus gigantes 