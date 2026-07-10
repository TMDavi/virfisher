import os

OUTDIR = config.get("outdir", "results")

rule fastp:
    input:
        forward=lambda wc: config["samples"][wc.sample]["forward"],
        reverseR=lambda wc: config["samples"][wc.sample]["reverseR"]
    output:
        forward=out("{sample}", "intermediate", "fastp", "{sample}_R1_trimmed.fastq.gz"),
        reverseR=out("{sample}", "intermediate","fastp", "{sample}_R2_trimmed.fastq.gz"),
        report=out("{sample}", "intermediate","fastp", "{sample}_fastp_report.html")

    log:
        stdout=out("{sample}", "intermediate","fastp", "stdout.log"),
        stderr=out("{sample}", "intermediate","fastp", "stderr.log")

    benchmark:
        out("{sample}", "intermediate", "fastp", "benchmark.txt")
    conda:
        "fastp"

    threads:
        config["resources"]["threads"]

    shell:
        """
        fastp \
            -i {input.forward} \
            -I {input.reverseR} \
            -o {output.forward} \
            -O {output.reverseR} \
            -q 20 \
            -w {threads} \
            -h {output.report} \
            1>{log.stdout} \
            2>{log.stderr}
        """