configfile: "config.yaml"
import os

rule all:
    input:
        #expand("results/{sample}/merged_viral_contigs.fasta", sample=config["samples"].keys())
        expand("results/{sample}/final_filtered/filtered_genomad.fasta", sample=config["samples"].keys())

rule virsorter2:
    input:
        assembly = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample])
    output:
        "results/{sample}/virsorter2/final-viral-combined.fa"
    params:
        outdir = "results/{sample}/virsorter2/"
    log:
        stdout = "results/{sample}/virsorter2/log-stdout.txt",
        stderr = "results/{sample}/virsorter2/log-stderr.txt"
    conda:
         "/MP_Data/mambaforge/envs/viral-id-sop"
    threads: 40
    shell:
        """
        virsorter run --keep-original-seq -i {input.assembly} -w {params.outdir} --include-groups dsDNAphage,ssDNA --min-length 5000 --min-score 0.5 -j {threads} all
        """

rule deepvirfinder:
    input:
        assembly = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample])
    output:
        "results/{sample}/deepvirfinder/{sample}_scaffolds.fasta_gt3000bp_dvfpred.txt"
    params:
        outdir = "results/{sample}/deepvirfinder/"
    log:
        stdout = "results/{sample}/deepvirfinder/log-stdout.txt",
        stderr = "results/{sample}/deepvirfinder/log-stderr.txt"
    conda:
         "/MP_Data/mambaforge/envs/dvf"
    threads: 50
    shell:
        """
        python /MP_Data/DeepVirFinder/dvf.py -i {input.assembly} -o {params.outdir} -l 3000 -c {threads}
        """

rule filter_dvf_contigs:
    input:
        assembly = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]),
        pred = "results/{sample}/deepvirfinder/{sample}_scaffolds.fasta_gt3000bp_dvfpred.txt"
    output:
        "results/{sample}/deepvirfinder/filtered_scaffolds.fasta"
    params:
        tempdir = lambda wildcards: f"results/{wildcards.sample}/deepvirfinder"
    shell:
        """
        python scripts/filter_dvf_contigs.py {input.assembly} {input.pred} {output} {params.tempdir}
        """

rule mergeviralcontigs:
    input:
        dvf = "results/{sample}/deepvirfinder/filtered_scaffolds.fasta",
        virsorter = "results/{sample}/virsorter2/final-viral-combined.fa"
    output:
        "results/{sample}/merged_viral/merged_viral_contigs.fasta"
    params:
        outdir = "results/{sample}/merged_viral_step1"
    shell:
        """
        python scripts/mergecontigs.py -o {params.outdir} -f {input.dvf} {input.virsorter}
        """

rule checkv:
    input:
        assembly = "results/{sample}/merged_viral_step1/merged_viral_contigs.fasta"
    output:
        "results/{sample}/checkv/viruses.fna",
        "results/{sample}/checkv/proviruses.fna",
        "results/{sample}/checkv/quality_summary.tsv"
    params:
        db = "/MP_Data/database/checkv-db-v1.5",
        outdir = "results/{sample}/checkv/"
    conda:
         "viral-id-sop"
    threads: 40
    shell:
        """
        checkv end_to_end {input} {params.outdir} -t {threads} -d {params.db}
        """

rule combine_checkv:
    input:
        viruses = "results/{sample}/checkv/viruses.fna",
        proviruses = "results/{sample}/checkv/proviruses.fna"
    output:
        "results/{sample}/checkv/combined.fna"
    shell:
        """
        python scripts/combine_checkv.py --virus {input.viruses} --provirus {input.proviruses} --output {output}
        """

rule filtering_one:
    input:
        assembly = "results/{sample}/checkv/combined.fna",
        quality = "results/{sample}/checkv/quality_summary.tsv"
    output:
        "results/{sample}/checkv/filtered_checkv.fasta"
    params:
        tempdir = "results/{sample}/checkv/tmp"
    shell:
        """
        python scripts/filtering_one.py --fasta_file {input.assembly} --quality_file {input.quality} --tempdir {params.tempdir} --output_file {output}
        """
        
rule genomad:
    input:
        filtered = "results/{sample}/checkv/filtered_checkv.fasta"
    output:
        "results/{sample}/genomad/filtered_checkv_annotate/filtered_checkv_genes.tsv"
    params:
        outdir = "results/{sample}/genomad/",
        db = "/MP_Data/database/genomad_db" #adcionar no config 
    log:
        stdout = "results/{sample}/genomad/log-stdout.txt",
        stderr = "results/{sample}/genomad/log-stderr.txt"
    conda:
         "/MP_Data/mambaforge/envs/genomad_env"
    threads: 40
    shell:
        """
        genomad end-to-end --cleanup {input.filtered} {params.outdir} {params.db} --threads {threads}
        """

rule filtering_two:
    input:
        assembly = "results/{sample}/checkv/filtered_checkv.fasta",
        annotation = "results/{sample}/genomad/filtered_checkv_annotate/filtered_checkv_genes.tsv"
    output:
        "results/{sample}/final_filtered/filtered_genomad.fasta"
    shell:
        """
        python scripts/filtering_two.py --fasta_file {input.assembly} --annotation_file {input.annotation} --output_file {output}
        """

rule final_quality:
    input:
        assembly = "results/{sample}/final_filtered/filtered_genomad.fasta",
        qual = "results/{sample}/checkv/quality_summary.tsv"
    output:
        "results/{sample}/final_filtered/final_quality_summary.tsv"
    shell:
        """
        python scripts/get_final_quality.py --fasta {input.assembly} --qual {input.qual} --output {output}
        """
rule final_annot_summary:
    input:
        assembly = "results/{sample}/final_filtered/filtered_genomad.fasta",
        annot = "results/{sample}/genomad/filtered_checkv_annotate/filtered_checkv_genes.tsv"
    output:
        "results/{sample}/final_filtered/final_annot_summary.tsv"
    shell:
        """
        python scripts/get_final_annot.py --fasta {input.assembly} --annot {input.annot} --output {output}
        """