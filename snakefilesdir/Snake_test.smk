configfile: "config.yaml"
import os

rule all:
    input:
        #expand("results/{sample}/final_filtered/final_quality_summary.tsv", sample=config["samples"].keys()),
        #expand("results/{sample}/final_filtered/final_annot_summary.tsv", sample=config["samples"].keys())
        expand("results/{sample}/final_filtered/with_hallmark/filtered_with_hallmark.fasta", sample=config["samples"].keys())


#rule final_quality:
#    input:
#        assembly = "results/{sample}/final_filtered/filtered_genomad.fasta",
#        qual = "results/{sample}/checkv/quality_summary.tsv"
#    output:
#        "results/{sample}/final_filtered/final_quality_summary.tsv"
#    shell:
#        """
#        python scripts/get_final_quality.py --fasta {input.assembly} --qual {input.qual} --output {output}
#        """


#rule final_annot_summary:
#    input:
#        assembly = "results/{sample}/final_filtered/filtered_genomad.fasta",
#        annot = "results/{sample}/genomad/filtered_checkv_annotate/filtered_checkv_genes.tsv"
#    output:
#        "results/{sample}/final_filtered/final_annot_summary.tsv"
#    shell:
#        """
#        python scripts/get_final_annot.py --fasta {input.assembly} --annot {input.annot} --output {output}
#        """

rule filter_hallmark:
    input:
        assembly = "results/{sample}/final_filtered/filtered_genomad.fasta",
        quality = "results/{sample}/final_filtered/final_quality_summary.tsv",
        annot = "results/{sample}/final_filtered/final_annot_summary.tsv"
    output:
        "results/{sample}/final_filtered/with_hallmark/filtered_with_hallmark.fasta",
        "results/{sample}/final_filtered/with_hallmark/filtered_with_hallmark_annotations.tsv",
        "results/{sample}/final_filtered/with_hallmark/filtered_with_hallmark_quality.tsv"
    params:
        outdir = "results/{sample}/final_filtered/with_hallmark"
    shell:
        """
        python scripts/filter_hallmark_phage.py --fasta_file {input.assembly} --annot_file {input.annot} --qual_file {input.quality} --out_dir {params.outdir}
        """