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