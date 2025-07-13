import os
configfile: 'config_read.yaml'



# Main target rule
rule all:
    input:
        "mapped_reads/idxstats_table.txt"

# Align reads to contigs using Bowtie2
rule bowtie:
    input:
        forward = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["forward"]),
        reverseR = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["reverseR"])
    output:
        temp("mapped_reads/{sample}/{sample}_vs_contigs.sam")
    log:
        "logs/bowtie2/{sample}.log"
    benchmark:
        "benchmarks/{sample}.bowtie2.benchmark.txt"
    threads: 50
    shell:
        "(bowtie2 -x bowtieDB/contigs -1 {input.forward} -2 {input.reverseR} "
        "-S {output} --local --no-unal -p {threads}) 2> {log}"

# Convert SAM to BAM
rule samtools_view:
    input:
        "mapped_reads/{sample}/{sample}_vs_contigs.sam"
    output:
        temp("mapped_reads/{sample}/{sample}_vs_contigs.bam")
    shell:
        "samtools view -Sb {input} >> {output}"

# Sort BAM
rule samtools_sort:
    input:
        "mapped_reads/{sample}/{sample}_vs_contigs.bam"
    output:
        temp("mapped_reads/{sample}/{sample}_vs_contigs_sorted.bam")
    shell:
        "samtools sort -o {output} {input}"

# Index sorted BAM (optional, temporary)
rule samtools_index:
    input:
        "mapped_reads/{sample}/{sample}_vs_contigs_sorted.bam"
    output:
        temp("mapped_reads/{sample}/{sample}_vs_contigs_sorted.bam.bai")
    shell:
        "samtools index {input}"

# Get idxstats from BAM
rule samtools_idxstats:
    input:
        "mapped_reads/{sample}/{sample}_vs_contigs_sorted.bam"
    output:
        "mapped_reads/{sample}/{sample}_vs_contigs_sorted.idxstats"
    shell:
        "samtools idxstats {input} >> {output}"

# Calculate per-contig coverage using CoverM
rule genome_coverage:
    input:
        "mapped_reads/{sample}/{sample}_vs_contigs_sorted.bam"
    output:
        "mapped_reads/{sample}/{sample}_contig_coverage.txt"
    conda:
        "/MP_Data/mambaforge/envs/coverm"
    shell:
        """
        coverm contig -b {input} -o {output}
        """
#JOin idxstats tables

rule join_idxstats:
    input:
        expand("mapped_reads/{sample}/{sample}_vs_contigs_sorted.idxstats",
               sample=config["samples"].keys())
    output:
        "mapped_reads/idxstats_table.txt"
    shell:
        """
        python scripts/process_map_data.py -i {input} -o {output}
        """
