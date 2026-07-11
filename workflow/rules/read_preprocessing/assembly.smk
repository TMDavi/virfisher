import os

OUTDIR = config.get("outdir", "results")

rule metaspades:
    input:
        forward = out("{sample}", "intermediate", "fastp", "{sample}_R1_trimmed.fastq.gz"),
        reverseR = out("{sample}", "intermediate", "fastp", "{sample}_R2_trimmed.fastq.gz")
    output:
       scaffolds = out("{sample}", "intermediate","metaspades", "scaffolds.fasta")
    params:
        outdir=out("{sample}", "intermediate","metaspades"),
        klist = config["metaspades"]["klist"]
    log:
        stdout = out("{sample}", "intermediate","metaspades", "stdout.log"),
        stderr = out("{sample}", "intermediate","metaspades", "stderr.log")
    conda:
         "metagenome"
    resources:
        mem_gb = config["metaspades"]["mem_gb"]
    threads:
        config["resources"]["threads"]
    shell:
        """
        metaspades.py -1 {input.forward} -2 {input.reverseR} \
        -o {params.outdir} -k {params.klist} \
        --threads {threads} --memory {resources.mem_gb} \
        > {log.stdout} \
        2> {log.stderr}
        """