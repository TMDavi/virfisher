import os

OUTDIR = config.get("outdir", "results")

def out(*parts):
    return os.path.join(OUTDIR, *parts)



include: "rules/qc.smk"
include: "rules/assembly.smk"

rule all:
    input:
        expand(out("{sample}", "fastp", "{sample}_fastp_report.html"), sample=config["samples"]),
        expand(out("{sample}", "metaspades", "scaffolds.fasta"), sample=config["samples"].keys())