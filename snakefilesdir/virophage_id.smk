configfile: 'config.yaml'

rule all:
    input:
        expand("results/{sample}/virophage_id/virophage_hmmer_results.domout", sample=config["samples"].keys()),
        expand("results/{sample}/virophage_id/virophage_parsed_results.txt", sample=config["samples"].keys()),
        expand("results/{sample}/virophage_id/virophage_filtered_results.txt", sample=config["samples"].keys()),
        expand("results/{sample}/virophage_id/virophage_extracted_genes.faa", sample=config["samples"].keys()),
        expand("results/{sample}/virophage_id/virophage_filtered_contigs.txt", sample=config["samples"].keys()),
        expand("results/{sample}/virophage_id/plv_hmmer_results.domout", sample=config["samples"].keys()),
        expand("results/{sample}/virophage_id/plv_parsed_results.txt", sample=config["samples"].keys()),
        expand("results/{sample}/virophage_id/plv_filtered_results.txt", sample=config["samples"].keys()),
        expand("results/{sample}/virophage_id/plv_extracted_genes.faa", sample=config["samples"].keys()),
        expand("results/{sample}/virophage_id/plv_filtered_contigs.txt", sample=config["samples"].keys())

rule virophage_extract:
    input:
        "results/{sample}/prodigal/{sample}_proteins.faa"
    output:
        "results/{sample}/virophage_id/virophage_hmmer_results.domout",
        "results/{sample}/virophage_id/virophage_parsed_results.txt",
        "results/{sample}/virophage_id/virophage_filtered_results.txt",
        "results/{sample}/virophage_id/virophage_extracted_genes.faa",
        "results/{sample}/virophage_id/virophage_filtered_contigs.txt",
        "results/{sample}/virophage_id/plv_hmmer_results.domout",
        "results/{sample}/virophage_id/plv_parsed_results.txt",
        "results/{sample}/virophage_id/plv_filtered_results.txt",
        "results/{sample}/virophage_id/plv_extracted_genes.faa",
        "results/{sample}/virophage_id/plv_filtered_contigs.txt"
    params:
        outdir = "results/{sample}/virophage_id"
    conda:
        "/MP_Data/mambaforge/envs/ncdlv_msearch"
    threads: 40
    shell:
        """
        python scripts/virophage_extract.py -i {input} -t {threads} -o {params.outdir} 
        """
rule filter_virophage_contigs:
    input:
        header_file = "results/{sample}/virophage_id/virophage_filtered_contigs.txt", 
        fasta_file = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample])
    params:
        temp_dir_parent = "results/{sample}/virophage_id"
    output:
        "results/{sample}/virophage_id/putative_virophage_contigs.fasta"
    threads: 20
    shell:
        """
        python scripts/filter_contigs.py -i {input.fasta_file} -hd {input.header_file} -o {output} --temp_dir_parent {params.temp_dir_parent} --num_chunks 20 -t {threads}
        """
rule filter_plv_contigs:
    input:
        header_file = "results/{sample}/virophage_id/plv_filtered_contigs.txt", 
        fasta_file = lambda wildcards: os.path.abspath(config["samples"][wildcards.sample])
    params:
        temp_dir_parent = "results/{sample}/virophage_id"
    output:
        "results/{sample}/virophage_id/putative_plv_contigs.fasta"
    threads: 20
    shell:
        """
        python scripts/filter_contigs.py -i {input.fasta_file} -hd {input.header_file} -o {output} --temp_dir_parent {params.temp_dir_parent} --num_chunks 20 -t {threads}
        """