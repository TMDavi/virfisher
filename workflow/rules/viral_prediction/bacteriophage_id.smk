rule virsorter2:
    input:
        scaffolds = out("{sample}", "intermediate","metaspades", "scaffolds.fasta")
    output:
        out("{sample}", "intermediate", "virsorter2", "final-viral-score.tsv")
    params:
        outdir = out("{sample}", "intermediate", "virsorter2")
    log:
        stdout = out("{sample}", "intermediate", "virsorter2", "log-stdout.txt"),
        stderr = out("{sample}", "intermediate", "virsorter2", "log-stderr.txt")
    conda:
         "viral-id-sop"
    threads: config["resources"]["threads"]
    shell:
        """
        virsorter run --keep-original-seq -i {input.scaffolds} -w {params.outdir} --include-groups dsDNAphage,ssDNA --min-length 1000 --min-score 0.5 -j {threads} all
        """
rule filter_virsorter2_contigs:
    input:
        score=out("{sample}", "intermediate", "virsorter2", "final-viral-score.tsv")
    output:
        ids=out("{sample}", "intermediate", "virsorter2", "viral_scaffolds.txt")
    shell:
        """
        tail -n +2 {input.score} | cut -f1 | sed 's/||.*//' > {output.ids}
        """

rule deepvirfinder:
    input:
        scaffolds = out("{sample}", "intermediate","metaspades", "scaffolds.fasta")
    output:
        out("{sample}", "intermediate","deepvirfinder","scaffolds.fasta_gt1000bp_dvfpred.txt")
    params:
        outdir = out("{sample}", "intermediate","deepvirfinder")
    log:
        stdout = out("{sample}", "intermediate", "deepvirfinder", "log-stdout.txt"),
        stderr = out("{sample}", "intermediate", "deepvirfinder", "log-stderr.txt")
    conda:
         "dvf"
    threads: config["resources"]["threads"]
    shell:
        """
        python programs/DeepVirFinder/dvf.py -i {input.scaffolds} -o {params.outdir} -l 1000 -c {threads}
        """

rule filter_dvf_contigs:
    input:
        pred=out("{sample}", "intermediate","deepvirfinder","scaffolds.fasta_gt1000bp_dvfpred.txt")
    output:
        ids=out("{sample}", "intermediate", "deepvirfinder", "viral_scaffolds.txt")
    shell:
        """
        tail -n +2 {input.pred} | cut -f1 > {output.ids}
        """

rule cenotetaker3:
    input:
        scaffolds=out("{sample}", "intermediate", "metaspades", "scaffolds.fasta")
    output:
        summary=out("{sample}","intermediate","cenote-taker3","{sample}","{sample}_virus_summary.tsv")
    params:
        outdir = out("{sample}", "intermediate", "cenote-taker3"),
        samplename = "{sample}"
    #log:
    #    stdout=out("{sample}", "intermediate", "cenote-taker3", "log-stdout.txt"),
    #    stderr=out("{sample}", "intermediate", "cenote-taker3", "log-stderr.txt")
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
        summary=out("{sample}","intermediate","cenote-taker3","{sample}","{sample}_virus_summary.tsv")
    output:
        ids=out("{sample}", "intermediate", "cenote-taker3", "viral_scaffolds.txt")
    shell:
        """
        tail -n +2 {input.summary} | cut -f2 > {output.ids}
        """

rule extract_viral_scaffolds:
    input:
        fasta=out("{sample}", "intermediate", "metaspades", "scaffolds.fasta"),
        vs2_ids=out("{sample}", "intermediate", "virsorter2", "viral_scaffolds.txt"),
        dvf_ids=out("{sample}", "intermediate", "deepvirfinder", "viral_scaffolds.txt"),
        ct3_ids=out("{sample}", "intermediate", "cenote-taker3", "viral_scaffolds.txt")
    output:
        fasta=out("{sample}", "intermediate", "viral_predicted_scaffolds_first_step.fasta"),
        ids=out("{sample}", "intermediate", "merged_ids.txt")
    shell:
        """
        cat {input.vs2_ids} {input.dvf_ids} {input.ct3_ids} | sort -u > {output.ids}

        seqkit grep -f {output.ids} {input.fasta} > {output.fasta}
        """

rule checkv:
    input:
        scaffolds = out("{sample}", "intermediate", "viral_predicted_scaffolds_first_step.fasta")
    output:
        viruses = out("{sample}","intermediate", "checkv", "viruses.fna"),
        proviruses = out("{sample}","intermediate", "checkv", "proviruses.fna"),
        summary = out("{sample}","intermediate", "checkv", "quality_summary.tsv")
    params:
        db = "databases/checkv-db-v1.5/", 
        outdir = out("{sample}","intermediate", "checkv")
    conda:
         "viral-id-sop"
    threads: config["resources"]["threads"]
    shell:
        """
        checkv end_to_end {input} {params.outdir} -t {threads} -d {params.db}
        """

rule combine_checkv:
    input:
        viruses = out("{sample}","intermediate", "checkv", "viruses.fna"),
        proviruses = out("{sample}","intermediate", "checkv", "proviruses.fna")
    output:
        out("{sample}", "intermediate", "checkv","combined.fna")
    shell:
        """
        python {WORKDIR}/scripts/combine_checkv.py --virus {input.viruses} --provirus {input.proviruses} --output {output}
        """

rule filtering_one:
    input:
        assembly = out("{sample}", "intermediate", "checkv","combined.fna"),
        quality = out("{sample}","intermediate", "checkv", "quality_summary.tsv")
    output:
        out("{sample}","intermediate","checkv","filtered_checkv.fasta")
    params:
        tempdir = out("{sample}","intermediate","checkv","tmp")
    shell:
        """
        python {WORKDIR}/scripts/filtering_one.py --fasta_file {input.assembly} --quality_file {input.quality} --tempdir {params.tempdir} --output_file {output}
        """
rule genomad:
    input:
        filtered = out("{sample}","intermediate","checkv","filtered_checkv.fasta")
    output:
        out("{sample}","intermediate", "genomad", "filtered_checkv_annotate", "filtered_checkv_genes.tsv")
    params:
        outdir = out("{sample}","intermediate", "genomad"),
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
        assembly = out("{sample}","intermediate","checkv","filtered_checkv.fasta"),
        annotation = out("{sample}","intermediate", "genomad", "filtered_checkv_annotate", "filtered_checkv_genes.tsv")
    output:
        out("{sample}","intermediate", "viral_predicted_scaffolds_second_step.fasta")
    shell:
        """
        python {WORKDIR}/scripts/filtering_two.py --fasta_file {input.assembly} --annotation_file {input.annotation} --output_file {output}
        """

