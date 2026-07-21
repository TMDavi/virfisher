import os

OUTDIR = config.get("outdir", "results")

def out(*parts):
    return os.path.join(OUTDIR, *parts)
    
include: "rules/read_preprocessing/qc.smk"
include: "rules/read_preprocessing/assembly.smk"
include: "rules/read_preprocessing/mapping_index.smk"
include: "rules/viral_prediction/bacteriophage_id.smk"
include: "rules/viral_prediction/NCLDV_virophage_plv_id.smk"


rule all:
    input:
        #QC
        expand(out("{sample}", "intermediate", "fastp", "{sample}_fastp_report.html"), sample=config["samples"]),

        #Assembly
        expand(out("{sample}", "intermediate", "metaspades", "scaffolds.fasta"), sample=config["samples"].keys()),

        #Bacteriophage prediction
        expand(out("{sample}", "intermediate", "viral_predicted_scaffolds_second_step.fasta"), sample=config["samples"].keys()),

        #Read mapping
        expand(
            out("{assembly}", "intermediate", "read_alignment",
                "{reads}_vs_{assembly}_sorted.bam"),
            assembly=config["samples"],
            reads=config["samples"]
        ),
        expand(out("{sample}", "final_results", "read_alignment", "{sample}_metabat.coverage"), sample=config["samples"].keys()),
        expand(out("{sample}", "final_results", "read_alignment", "{sample}.coverage"), sample=config["samples"].keys()),

        #NCLDV identification
        expand(out("{sample}", "intermediate","BEREN","Final_results","Run_Summary.txt"), sample=config["samples"].keys()),

        #expand(out("{sample}", "intermediate","prodigal","{sample}_proteins.faa"), sample=config["samples"].keys()),
        #expand(out("{sample}","intermediate","virophage_plv_id","genomad","viruses_annotate","viruses_genes.tsv"), sample=config["samples"].keys())
