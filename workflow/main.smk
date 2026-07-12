import os

OUTDIR = config.get("outdir", "results")

def out(*parts):
    return os.path.join(OUTDIR, *parts)



include: "rules/read_preprocessing/qc.smk"
include: "rules/read_preprocessing/assembly.smk"
include: "rules/read_preprocessing/alignment_DB_bwa.smk"
include: "rules/read_preprocessing/read_alignment.smk"
include: "rules/viral_prediction/bacteriophage_id.smk"

rule all:
    input:
        expand(out("{sample}", "intermediate", "fastp", "{sample}_fastp_report.html"), sample=config["samples"]),
        expand(out("{sample}", "intermediate", "metaspades", "scaffolds.fasta"), sample=config["samples"].keys()),
        expand(out("mappingDB", "mappingDB.{suffix}"),suffix=["amb","ann","bwt","pac","sa"]),
        out("final_results", "read_mapping", "merged_idxstats.tsv"),
        out("final_results", "read_mapping", "merged_depth.tsv"),
        out("final_results", "read_mapping", "merged_coverage.tsv"),
        expand(out("{sample}", "intermediate", "viral_predicted_scaffolds_second_step.fasta"), sample=config["samples"].keys())
