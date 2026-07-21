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
        out("{sample}", "intermediate","prodigal","{sample}_proteins.faa")
    output:
        out("{sample}", "intermediate", "ncldv_id", "hmmer_results.domout"),
        out("{sample}", "intermediate", "ncldv_id", "parsed_results.txt"),
        out("{sample}", "intermediate", "ncldv_id", "filtered_results.txt"),
        out("{sample}", "intermediate", "ncldv_id", "extracted_genes.faa"),
        out("{sample}", "intermediate", "ncldv_id", "putative_ncldv_contig_list.txt")
    params:
        outdir = out("{sample}", "intermediate", "ncldv_id")
    conda:
        "ncdlv_msearch"
    threads: config["resources"]["threads"]
    shell:
        "python {WORKDIR}/scripts/ncldv_extract.py -i {input} -t {threads} -o {params.outdir}"

rule extract_viral_scaffolds:
    input:
        header_file = "results/{sample}/ncldv_id/putative_ncldv_contig_list.txt", 
        fasta_file = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample])
    output:
        fasta=out("{sample}", "intermediate", "viral_predicted_scaffolds_first_step.fasta"),
        ids=out("{sample}", "intermediate", "merged_ids.txt")
    shell:
        """
        cat {input.vs2_ids} {input.dvf_ids} {input.ct3_ids} | sort -u > {output.ids}

        seqkit grep -f {output.ids} {input.fasta} > {output.fasta}
        """

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