configfile: "config.yaml"
import os

rule all:
    input:
        expand("results/{sample}/ncldv_id/hmmer_results.domout", sample=config['samples'].keys()),
        expand("results/{sample}/ncldv_id/parsed_results.txt", sample=config['samples'].keys()),
        expand("results/{sample}/ncldv_id/filtered_results.txt", sample=config['samples'].keys()),
        expand("results/{sample}/ncldv_id/extracted_genes.faa", sample=config['samples'].keys()),
        expand("results/{sample}/ncldv_id/putative_ncldv_contig_list.txt", sample=config['samples'].keys()),
        expand("results/{sample}/ncldv_id/putative_ncldv_contigs.fasta", sample=config['samples'].keys()),
        expand("results/{sample}/ncldv_id/final_ncldv_contigs.fasta", sample=config['samples'].keys()),
        

rule ncdlv_extract:
    input:
        "results/{sample}/prodigal/{sample}_proteins.faa"
    output:
        "results/{sample}/ncldv_id/hmmer_results.domout",
        "results/{sample}/ncldv_id/parsed_results.txt",
        "results/{sample}/ncldv_id/filtered_results.txt",
        "results/{sample}/ncldv_id/extracted_genes.faa",
        "results/{sample}/ncldv_id/putative_ncldv_contig_list.txt"
    params:
        outdir = "results/{sample}/ncldv_id"
    conda:
        "/MP_Data/mambaforge/envs/ncdlv_msearch"
    threads: 40
    shell:
        "python scripts/ncldv_extract.py -i {input} -t {threads} -o {params.outdir}"

rule filter_ncldv_contigs:
    input:
        header_file = "results/{sample}/ncldv_id/putative_ncldv_contig_list.txt", 
        fasta_file = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample])
    params:
        temp_dir_parent = "results/{sample}/ncldv_id"
    output:
        "results/{sample}/ncldv_id/putative_ncldv_contigs.fasta"
    threads: 20
    shell:
        """
        python scripts/filter_contigs.py -i {input.fasta_file} -hd {input.header_file} -o {output} --temp_dir_parent {params.temp_dir_parent} --num_chunks 20 -t {threads}
        """
rule viralrecall:
    input:
        "results/{sample}/ncldv_id/putative_ncldv_contigs.fasta"
    output:
        "results/{sample}/ncldv_id/viralrecall/{sample}.summary.tsv"
    params:
        samplename = lambda wildcards: wildcards.sample,
        outdir = "results/{sample}/ncldv_id/viralrecall"
    conda:
        "/MP_Data/mambaforge/envs/viralrecall"
    threads: 40
    shell:
        """
        python programs/viralrecall/viralrecall.py -i {input} -p {params.samplename} -t {threads} -c -f -o {params.outdir}
        """

rule get_contigs_positive_score:
    input:
        "results/{sample}/ncldv_id/viralrecall/{sample}.summary.tsv"
    output:
        "results/{sample}/ncldv_id/viralrecall/ncldv_positive_score.txt"
    shell:
        """
        python scripts/get_pos_score_ncldv.py -i {input} -o {output}
        """

rule filter_ncldv_contigs2:
    input:
        header_file = "results/{sample}/ncldv_id/viralrecall/ncldv_positive_score.txt", 
        fasta_file = "results/{sample}/ncldv_id/putative_ncldv_contigs.fasta"
    params:
        temp_dir_parent = "results/{sample}/ncldv_id"
    output:
        "results/{sample}/ncldv_id/final_ncldv_contigs.fasta"
    threads: 20
    shell:
        """
        python scripts/filter_contigs.py -i {input.fasta_file} -hd {input.header_file} -o {output} --temp_dir_parent {params.temp_dir_parent} --num_chunks 20 -t {threads}
        """