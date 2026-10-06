import os

OUTDIR = config.get("outdir", "results")

def out(*parts):
    return os.path.join(OUTDIR, *parts)
    
include: "rules/read_preprocessing/qc.smk"
include: "rules/read_preprocessing/assembly.smk"
include: "rules/read_preprocessing/read_alignment.smk"
include: "rules/viral_prediction/bacteriophage_id.smk"
include: "rules/dereplication/dereplication.smk"
include: "rules/host_prediction_and_lifestyle/iphop.smk"
include: "rules/host_prediction_and_lifestyle/vibrant.smk"
include: "rules/host_prediction_and_lifestyle/bacphlip.smk"
include: "rules/quality/checkv.smk"
include: "rules/annotation/pharokka.smk"
include: "rules/taxonomy/genomad.smk"
include: "rules/taxonomy/vcontact3.smk"
include: "rules/annotation/dramv.smk"

rule all:
    input:
        expand(out("per_sample_results", "{sample}", "intermediate", "fastp", "{sample}_fastp_report.html"),
            sample=config["samples"].keys()),

        expand(out("per_sample_results", "{sample}", "intermediate", "metaspades", "scaffolds.fasta"), 
            sample=config["samples"].keys()),

        #Bacteriophage prediction
        expand(out("per_sample_results", "{sample}","intermediate", "final_phage_sequences.fasta"), 
            sample=config["samples"].keys()),

        #Read mapping
        expand(
            out("per_sample_results","{assembly}", "intermediate", "read_alignment",
                "{reads}_vs_{assembly}_sorted.bam"),
            assembly=config["samples"],
            reads=config["samples"]
        ),
        expand(out("per_sample_results", "{sample}", "final_results", "read_alignment", "{sample}_metabat.coverage"), sample=config["samples"].keys()),
        expand(out("per_sample_results", "{sample}", "final_results", "read_alignment", "{sample}.coverage"), sample=config["samples"].keys()),

        #Multiqc

        #Dereplication
        expand(out("per_sample_results","{sample}","final_results","dereplicated","miuvigs.fasta"), sample=config['samples'].keys()),

        #Host predicition
        #Iphop
        expand(out("per_sample_results", "{sample}","final_results","Host_prediction","Host_prediction_to_genome_m90.csv"), sample=config['samples'].keys()),

        #Life style 
        expand(out("per_sample_results", "{sample}","final_results","life_style","vibrant_results", "final_phage_sequences.phages_lysogenic.fna"),sample=config["samples"].keys()),
        expand(out("per_sample_results", "{sample}","final_results","life_style","vibrant_results", "final_phage_sequences.phages_lytic.fna"),sample=config["samples"].keys()),
        expand(out("per_sample_results", "{sample}","final_results","life_style","bacphlip"),sample=config["samples"].keys()),

        #Annotation
        #DRAMv

        #pharaokka
        expand(out("per_sample_results", "{sample}","final_results", "annotation","pharokka","pharokka_cds_final_merged_output.tsv"),
            sample=config["samples"].keys()),

        #Miuvig Quality
        #Checkv
        expand(out("per_sample_results", "{sample}","final_results", "quality_summary", "quality_summary.tsv"),
            sample=config["samples"].keys()),

        #Taxonomy
        #Genomad
        expand(out("per_sample_results", "{sample}","final_results", "taxonomy","genomad", "final_phage_sequences_summary", "final_phage_sequences_virus_summary.tsv"),
            sample=config["samples"].keys()),
        #Vcontact3
        out("merged_results","vcontact3","exports","final_assignments.csv")

        