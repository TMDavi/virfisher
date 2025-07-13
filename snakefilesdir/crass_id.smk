configfile: "config.yaml"
import os

rule all:
    input:
        expand("results/{sample}/crass_id/hmmer_results.domout", sample=config["samples"].keys()),
        expand("results/{sample}/crass_id/parsed_results.txt", sample=config["samples"].keys()),
        expand("results/{sample}/crass_id/filtered_results.txt", sample=config["samples"].keys()),
        expand("results/{sample}/crass_id/extracted_genes.faa", sample=config["samples"].keys()),
        expand("results/{sample}/crass_id/putative_crass_contig_list.txt", sample=config["samples"].keys()),
        expand("results/{sample}/crass_id/putative_crass_contigs.fasta", sample=config["samples"].keys())

rule crass_extract:
    input:
        "results/{sample}/prodigal/{sample}_proteins.faa"
    output:
        "results/{sample}/crass_id/hmmer_results.domout",
        "results/{sample}/crass_id/parsed_results.txt",
        "results/{sample}/crass_id/filtered_results.txt",
        "results/{sample}/crass_id/extracted_genes.faa",
        "results/{sample}/crass_id/putative_crass_contig_list.txt"
    params:
        outdir = "results/{sample}/crass_id/"
    conda:
        "/MP_Data/mambaforge/envs/ncdlv_msearch"
    threads: 40
    shell:
        "python scripts/crass_extract.py -i {input} -t {threads} -o {params.outdir}"

rule filter_crass_contigs:
    input:
        header_file = "results/{sample}/crass_id/putative_crass_contig_list.txt", 
        fasta_file = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample])
    params:
        temp_dir_parent = "results/{sample}/crass_id/"
    output:
        "results/{sample}/crass_id/putative_crass_contigs.fasta"
    threads: 20
    shell:
        """
        python scripts/filter_contigs.py -i {input.fasta_file} -hd {input.header_file} -o {output} --temp_dir_parent {params.temp_dir_parent} --num_chunks 20 -t {threads}
        """