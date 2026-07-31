rule virsorter2:
    input:
        scaffolds = out("per_sample_results", "{sample}", "intermediate","metaspades", "scaffolds.fasta")
    output:
        out("per_sample_results", "{sample}", "intermediate", "virsorter2", "final-viral-score.tsv")
    params:
        outdir = out("per_sample_results", "{sample}", "intermediate", "virsorter2")
    log:
        stdout = out("per_sample_results", "{sample}", "intermediate", "virsorter2", "log-stdout.txt"),
        stderr = out("per_sample_results", "{sample}", "intermediate", "virsorter2", "log-stderr.txt")
    conda:
         "viral-id-sop"
    threads: config["resources"]["threads"]
    shell:
        """
        virsorter run --keep-original-seq -i {input.scaffolds} -w {params.outdir} --include-groups dsDNAphage,ssDNA --min-length 1000 --min-score 0.5 -j {threads} all
        """
rule filter_virsorter2_contigs:
    input:
        score=out("per_sample_results", "{sample}", "intermediate", "virsorter2", "final-viral-score.tsv")
    output:
        ids=out("per_sample_results", "{sample}", "intermediate", "virsorter2", "viral_scaffolds.txt")
    shell:
        """
        tail -n +2 {input.score} | cut -f1 | sed 's/||.*//' > {output.ids}
        """

rule deepvirfinder:
    input:
        scaffolds = out("per_sample_results", "{sample}", "intermediate","metaspades", "scaffolds.fasta")
    output:
        out("per_sample_results", "{sample}", "intermediate","deepvirfinder","scaffolds.fasta_gt1000bp_dvfpred.txt")
    params:
        outdir = out("per_sample_results", "{sample}", "intermediate","deepvirfinder")
    log:
        stdout = out("per_sample_results", "{sample}", "intermediate", "deepvirfinder", "log-stdout.txt"),
        stderr = out("per_sample_results", "{sample}", "intermediate", "deepvirfinder", "log-stderr.txt")
    conda:
         "dvf"
    threads: config["resources"]["threads"]
    shell:
        """
        python programs/DeepVirFinder/dvf.py -i {input.scaffolds} -o {params.outdir} -l 1000 -c {threads}
        """

rule filter_dvf_contigs:
    input:
        pred=out("per_sample_results", "{sample}", "intermediate","deepvirfinder","scaffolds.fasta_gt1000bp_dvfpred.txt")
    output:
        ids=out("per_sample_results", "{sample}", "intermediate", "deepvirfinder", "viral_scaffolds.txt")
    shell:
        """
        tail -n +2 {input.pred} | cut -f1 > {output.ids}
        """

rule cenotetaker3:
    input:
        scaffolds=out("per_sample_results", "{sample}", "intermediate", "metaspades", "scaffolds.fasta")
    output:
        summary=out("per_sample_results", "{sample}","intermediate","cenote-taker3","{sample}","{sample}_virus_summary.tsv")
    params:
        outdir = out("per_sample_results", "{sample}", "intermediate", "cenote-taker3"),
        samplename = "{sample}"
    #log:
    #    stdout=out("per_sample_results", "{sample}", "intermediate", "cenote-taker3", "log-stdout.txt"),
    #    stderr=out("per_sample_results", "{sample}", "intermediate", "cenote-taker3", "log-stderr.txt")
    conda:
         "ct3_env"
    threads:
        config["resources"]["threads"]
    shell:
        """
        cenotetaker3 -c {input.scaffolds} -r {params.samplename} -wd {params.outdir} -p T -db virion dnarep -t {threads}
        """

rule filter_cenote_taker_contigs:
    input:
        summary=out("per_sample_results", "{sample}","intermediate","cenote-taker3", "{sample}","{sample}_virus_summary.tsv")
    output:
        ids=out("per_sample_results", "{sample}", "intermediate", "cenote-taker3", "viral_scaffolds.txt")
    shell:
        """
        tail -n +2 {input.summary} | cut -f2 > {output.ids}
        """

rule extract_viral_scaffolds:
    input:
        fasta=out("per_sample_results", "{sample}", "intermediate", "metaspades", "scaffolds.fasta"),
        vs2_ids=out("per_sample_results", "{sample}", "intermediate", "virsorter2", "viral_scaffolds.txt"),
        dvf_ids=out("per_sample_results", "{sample}", "intermediate", "deepvirfinder", "viral_scaffolds.txt"),
        ct3_ids=out("per_sample_results", "{sample}", "intermediate", "cenote-taker3", "viral_scaffolds.txt")
    output:
        fasta=out("per_sample_results", "{sample}", "intermediate", "viral_predicted_scaffolds_first_step.fasta"),
        ids=out("per_sample_results", "{sample}", "intermediate", "merged_ids.txt")
    shell:
        """
        cat {input.vs2_ids} {input.dvf_ids} {input.ct3_ids} | sort -u > {output.ids}

        seqkit grep -f {output.ids} {input.fasta} > {output.fasta}
        """

rule checkv:
    input:
        scaffolds = out("per_sample_results", "{sample}", "intermediate", "viral_predicted_scaffolds_first_step.fasta")
    output:
        viruses = out("per_sample_results", "{sample}","intermediate", "checkv", "viruses.fna"),
        proviruses = out("per_sample_results", "{sample}","intermediate", "checkv", "proviruses.fna"),
        summary = out("per_sample_results", "{sample}","intermediate", "checkv", "quality_summary.tsv")
    params:
        db = "databases/checkv-db-v1.5/", 
        outdir = out("per_sample_results", "{sample}","intermediate", "checkv")
    conda:
         "viral-id-sop"
    threads: config["resources"]["threads"]
    shell:
        """
        checkv end_to_end {input} {params.outdir} -t {threads} -d {params.db}
        """

rule combine_checkv:
    input:
        viruses = out("per_sample_results", "{sample}","intermediate", "checkv", "viruses.fna"),
        proviruses = out("per_sample_results", "{sample}","intermediate", "checkv", "proviruses.fna")
    output:
        out("per_sample_results", "{sample}", "intermediate", "checkv","combined.fna")
    shell:
        """
        python {WORKDIR}/scripts/combine_checkv.py --virus {input.viruses} --provirus {input.proviruses} --output {output}
        """

rule filtering_one:
    input:
        assembly = out("per_sample_results", "{sample}", "intermediate", "checkv","combined.fna"),
        quality = out("per_sample_results", "{sample}","intermediate", "checkv", "quality_summary.tsv")
    output:
        out("per_sample_results", "{sample}","intermediate","checkv","filtered_checkv.fasta")
    params:
        tempdir = out("per_sample_results", "{sample}","intermediate","checkv","tmp")
    shell:
        """
        python {WORKDIR}/scripts/filtering_one.py --fasta_file {input.assembly} --quality_file {input.quality} --tempdir {params.tempdir} --output_file {output}
        """
rule genomad:
    input:
        filtered = out("per_sample_results", "{sample}","intermediate","checkv","filtered_checkv.fasta")
    output:
        out("per_sample_results", "{sample}","intermediate", "genomad", "filtered_checkv_annotate", "filtered_checkv_genes.tsv")
    params:
        outdir = out("per_sample_results", "{sample}","intermediate", "genomad"),
        db = "databases/genomad_db" 
    conda:
         "genomad_env"
    threads: config["resources"]["threads"]
    shell:
        """
        genomad end-to-end --cleanup {input.filtered} {params.outdir} {params.db} --threads {threads}
        """

rule filtering_two:
    input:
        assembly = out("per_sample_results", "{sample}","intermediate","checkv","filtered_checkv.fasta"),
        annotation = out("per_sample_results", "{sample}","intermediate", "genomad", "filtered_checkv_annotate", "filtered_checkv_genes.tsv")
    output:
        out("per_sample_results", "{sample}","intermediate", "final_phage_sequences.fasta")
    shell:
        """
        python {WORKDIR}/scripts/filtering_two.py --fasta_file {input.assembly} --annotation_file {input.annotation} --output_file {output}
        """

