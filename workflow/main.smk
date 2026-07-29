import os

OUTDIR = config.get("outdir", "results")

def out(*parts):
    return os.path.join(OUTDIR, *parts)
    
include: "rules/read_preprocessing/qc.smk"
include: "rules/read_preprocessing/assembly.smk"
include: "rules/read_preprocessing/read_alignment.smk"
include: "rules/viral_prediction/bacteriophage_id.smk"
include: "rules/dereplication/dereplication.smk"
#include: "rules/host_prediction_and_lifestyle/iphop.smk"
include: "rules/host_prediction_and_lifestyle/vibrant.smk"
include: "rules/host_prediction_and_lifestyle/bacphlip.smk"

rule all:
    input:
        #Setup databases

        #QC
        expand(out("{sample}", "intermediate", "fastp", "{sample}_fastp_report.html"), sample=config["samples"]),

        #Assembly
        expand(out("{sample}", "intermediate", "metaspades", "scaffolds.fasta"), sample=config["samples"].keys()),

        #Bacteriophage prediction
        expand(out("{sample}","intermediate", "final_phage_sequences.fasta"), sample=config["samples"].keys()),

        #Read mapping
        expand(
            out("{assembly}", "intermediate", "read_alignment",
                "{reads}_vs_{assembly}_sorted.bam"),
            assembly=config["samples"],
            reads=config["samples"]
        ),
        expand(out("{sample}", "final_results", "read_alignment", "{sample}_metabat.coverage"), sample=config["samples"].keys()),
        expand(out("{sample}", "final_results", "read_alignment", "{sample}.coverage"), sample=config["samples"].keys()),

        #Multiqc

        #Dereplication
        expand(out("{sample}","final_results","dereplicated","miuvigs.fasta"), sample=config['samples'].keys()),

        #Host predicition
        #Iphop

        #Life style 
        expand(out("{sample}","final_results","life_style","vibrant_results", "final_phage_sequences.phages_lysogenic.fna"),sample=config["samples"].keys()),
        expand(out("{sample}","final_results","life_style","vibrant_results", "final_phage_sequences.phages_lytic.fna"),sample=config["samples"].keys()),
        expand(out("{sample}","final_results","life_style","bacphlip"),sample=config["samples"].keys())

        #Annotation
        #DRAMv
        #pharaokka

        #Miuvig Quality
        #Checkv

        #Taxonomy
        #Genomad
        #Vcontact3

        