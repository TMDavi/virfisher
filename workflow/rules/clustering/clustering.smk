configfile: 'config.yaml'

rule merge_contigs:
    input:
        phages="results/{sample}/final_filtered/filtered_genomad.fasta"
        crass="results/{sample}/crass_id/putative_crass_contigs.fasta"
        ncldv="results/{sample}/ncldv_id/putative_ncldv_contigs.fasta"
        virophage="results/{sample}/virophage_id/putative_virophage_contigs.fasta"
        plv="results/{sample}/virophage_id/putative_plv_contigs.fasta"
    output:
        
    params:
        outdir = "results/{sample}/merged_viral_step2"