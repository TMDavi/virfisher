configfile: "config.yaml"
import os

rule all:
    input:
        expand("results/{sample}/fastp/{sample}_fastp_report.html", sample=config["samples"].keys()),
        expand("results/{sample}/assembly/metaspades/scaffolds.fasta", sample=config["samples"].keys()),
        expand("results/{sample}/merged_viral_contigs.fasta", sample=config["samples"].keys()),
        expand("results/{sample}/checkv/.complete", sample=config["samples"].keys()),
        expand("results/{sample}/checkv/combined.fna", sample=config["samples"].keys()),
        expand("results/{sample}/checkv/filtered_checkv.fasta", sample=config["samples"].keys()),
        expand("results/{sample}/genomad/contigs.fasta", sample=config["samples"].keys()),
        expand("results/{sample}/final_filtered/filtered_genomad.fasta", sample=config["samples"].keys())

rule fastp:
    input:
        forward = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["forward"]),
        reverseR = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample]["reverseR"])
    output:
        forward = "results/{sample}/fastp/{sample}_R1_trimmed.fastq.gz",
        reverseR = "results/{sample}/fastp/{sample}_R2_trimmed.fastq.gz",
        report = "results/{sample}/fastp/{sample}_fastp_report.html"
    params:
        outdir = "results/{sample}/fastp"
    log:
        stdout = "results/{sample}/fastp/log-stdout.txt",
        stderr = "results/{sample}/fastp/log-stderr.txt"
    conda:
        "envs/fastp.yaml"
    benchmark:
        "results/{sample}/fastqc/benchmark.txt"
    threads:
        config["threads"]
    shell:
        """
        fastp -i {input.forward} -I {input.reverseR} \ 
        -o {output.forward} -O {output.reverseR}  \ 
        -q 20 -w 15 \
        -h {output.report} 
        """
    
#rule megahit:
#    input:
#        forward = "results/{sample}/fastp/{sample}_R1_trimmed.fastq.gz",
#        reverseR = "results/{sample}/fastp/{sample}_R2_trimmed.fastq.gz"
#    params:
#        klist = config["megahit"]['klist'],
#        output = "results/{sample}/assembly/megahit"
#    output:
#        "results/{sample}/assembly/megahit/final.contigs.fa"
#    log:
#        stdout = "results/{sample}/assembly/megahit/log-stdout.txt",
#        stderr = "results/{sample}/assembly/megahit/log-stderr.txt"
#    conda:
#        "envs/megahit.yaml"
#    resources:
#        mem_gb = 100
#    threads:
#        config["threads"]
#    shell:
#        """
#            megahit -f -1 {input[0]} -2 {input[1]} -t {threads} --presets meta-large -o {params.output} --min-contig-len 300
#        """ 

rule metaspades:
    input:
        forward = "results/{sample}/fastp/{sample}_R1_trimmed.fastq.gz",
        reverseR = "results/{sample}/fastp/{sample}_R2_trimmed.fastq.gz"
    output:
        "results/{sample}/assembly/metaspades/scaffolds.fasta"
    params:
        klist = "27,37,47,57,67,77,87,97,107,117,127",
        outdir = "results/{sample}/assembly/metaspades"
    log:
        stdout = "results/{sample}/assembly/metaspades/log-stdout.txt",
        stderr = "results/{sample}/assembly/metaspades/log-stderr.txt"
    conda:
         "envs/spades.yaml"
    resources:
        mem_gb = 500
    threads:
        config["threads"]
    shell:
        """
        metaspades.py -1 {input.forward} -2 {input.reverseR} -o {params.outdir} -k {params.klist} --threads {threads} --memory {resources.mem_gb}
        """

rule virsorter2:
    input:
        assembly = "results/{sample}/assembly/metaspades/scaffolds.fasta"
    output:
        "results/{sample}/virsorter2/final-viral-combined.fa"
    params:
        outdir = "results/{sample}/virsorter2/"
    log:
        stdout = "results/{sample}/virsorter2/log-stdout.txt",
        stderr = "results/{sample}/virsorter2/log-stderr.txt"
    conda:
         "envs/virsorter2.yaml"
    threads:
        config["threads"]
    shell:
        """
        virsorter run --keep-original-seq -i {input.assembly} -w {params.outdir} --include-groups dsDNAphage,ssDNA --min-length 5000 --min-score 0.5 -j {threads} all
        """

rule deepvirfinder:
    input:
        assembly = "results/{sample}/assembly/metaspades/scaffolds.fasta"
    output:
        "results/{sample}/deepvirfinder/scaffolds.fasta_gt3000bp_dvf_pred.txt"
    params:
        outdir = "results/{sample}/deepvirfinder/"
    log:
        stdout = "results/{sample}/deepvirfinder/log-stdout.txt",
        stderr = "results/{sample}/deepvirfinder/log-stderr.txt"
    conda:
         "envs/deepvirfinder.yaml"
    threads:
        config["threads"]
    shell:
        """
        python dvf.py -i {params.assembly} -o {params.outdir} -l 3000 -c {threads}
        """

#rule cenotetaker3:
#    input:
#        assembly = "results/{sample}/assembly/metaspades/scaffolds.fasta"
#    output:
#        viral_seqs = "results/{sample}/cenote-taker3/{sample}/{sample}_virus_sequences.fna"
#    params:
#        outdir = "results/{sample}/cenote-taker3/",
#        samplename = "{sample}"
#    log:
#        stdout = "results/{sample}/cenote-taker3/log-stdout.txt",
#        stderr = "results/{sample}/cenote-taker3/log-stderr.txt"
#    conda:
#         "envs/cenote-taker3.yaml"
#    threads:
#        config["threads"]
#    shell:
#        """
#        cd {params.outdir}
#        cenotetaker3 -c {input.assembly} -r {params.samplename} -p T -db virion dnarep -t {threads}
#        cd {config["workdir"]}
#        """

rule filter_dvf_contigs:
    input:
        assembly = "results/{sample}/assembly/metaspades/scaffolds.fasta",
        pred = "results/{sample}/deepvirfinder/scaffolds.fasta_gt3000bp_dvf_pred.txt"
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
        "results/{sample}/merged_viral_contigs.fasta"
    shell:
        """
        python scripts/mergecontigs.py --output {output} --dvf_file {input.dvf} --vs2_file {input.virsorter}
        """
        
rule checkv:
    input:
        assembly = "results/{sample}/merged_viral_contigs.fasta"
    output:
        "results/{sample}/checkv/viruses.fna",
        "results/{sample}/checkv/proviruses.fna",
        "results/{sample}/checkv/quality_summary.tsv"
    params:
        db = "/MP_Data/database/checkv-db-v1.5",
        outdir = "results/{sample}/checkv/"
    conda:
         "envs/checkv.yaml"
    threads:
        config["threads"]
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
         "envs/genomad.yaml"
    threads:
        config["threads"]
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