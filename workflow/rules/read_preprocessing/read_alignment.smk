WORKDIR = config["WORKDIR"]

rule bwa_mem:
    input:
        index=expand(
            out("mappingDB", "mappingDB.{suffix}"),
            suffix=["amb", "ann", "bwt", "pac", "sa"],
        ),
        forward=out("{sample}", "intermediate", "fastp", "{sample}_R1_trimmed.fastq.gz"),
        reverseR=out("{sample}", "intermediate", "fastp", "{sample}_R2_trimmed.fastq.gz"),
    output:
        bam=out("{sample}", "intermediate", "mapped_reads", "{sample}_vs_" + "all_scaffolds" + "_sorted.bam")
    log:
        stderr=out("{sample}", "logs", "bwa.stderr.log")
    benchmark:
        out("{sample}", "benchmarks", "bwa.txt")
    threads:
        config["resources"]["threads"]
    params:
        index_prefix=out("mappingDB", "mappingDB")
    conda:
        "coverm"
    shell:
        """
        bwa mem -t {threads} {params.index_prefix} {input.forward} {input.reverseR} 2> {log.stderr} | samtools sort -@ {threads} -o {output.bam}
        """

rule samtools_index:
    input:
        bam=out("{sample}", "intermediate", "mapped_reads", "{sample}_vs_" + "all_scaffolds" + "_sorted.bam")
    output:
        bai=out("{sample}", "intermediate", "mapped_reads", "{sample}_vs_" + "all_scaffolds" + "_sorted.bam.bai")
    threads:
        config["resources"]["threads"]
    conda:
        "coverm"
    shell:
        """
        samtools index -@ {threads} {input.bam} {output.bai}
        """

rule samtools_idxstats:
    input:
        bam=out("{sample}", "intermediate", "mapped_reads", "{sample}_vs_" + "all_scaffolds" + "_sorted.bam")
    output:
        out("{sample}", "intermediate", "mapped_reads", "{sample}_vs_" + "all_scaffolds" + "_idxstats.tsv")
    conda:
        "coverm"
    shell:
        """
        samtools idxstats {input.bam} > {output}
        """

rule genome_depth:
    input:
        bam=out("{sample}", "intermediate", "mapped_reads", "{sample}_vs_" + "all_scaffolds" + "_sorted.bam")
    output:
        out("{sample}", "intermediate", "mapped_reads", "{sample}_vs_" + "all_scaffolds" + "_depth.tsv")
    conda:
        "coverm"
    shell:
        """
        coverm contig -b {input.bam} -o {output}
        """

rule calc_coverage:
    input:
        bam=out("{sample}", "intermediate", "mapped_reads", "{sample}_vs_" + "all_scaffolds" + "_sorted.bam")
    output:
        out("{sample}", "intermediate", "mapped_reads", "{sample}_vs_" + "all_scaffolds" + "_coverage.tsv")
    conda:
        "coverm"
    shell:
        """
        bedtools genomecov -ibam {input.bam} > {output}
        """

rule merge_idxstats:
    input:
        expand(
            out("{sample}", "intermediate", "mapped_reads", "{sample}_vs_" + "all_scaffolds" + "_idxstats.tsv"),
            sample=config["samples"]
        )
    output:
        out("final_results", "read_mapping", "merged_idxstats.tsv")
    conda:
        "coverm"
    shell:
        """
        python {WORKDIR}/scripts/merge_idxstats.py --files {input} --output {output}
        """

rule merge_depth:
    input:
        expand(
            out("{sample}", "intermediate", "mapped_reads", "{sample}_vs_" + "all_scaffolds" + "_depth.tsv"),
            sample=config["samples"]
        )
    output:
        out("final_results", "read_mapping", "merged_depth.tsv")
    conda:
        "coverm"
    shell:
        """
        python {WORKDIR}/scripts/merge_depth.py --files {input} --output {output}
        """

rule merge_coverage:
    input:
        expand(
            out("{sample}", "intermediate", "mapped_reads", "{sample}_vs_" + "all_scaffolds" + "_coverage.tsv"),
            sample=config["samples"]
        )
    output:
        out("final_results", "read_mapping", "merged_coverage.tsv")
    conda:
        "coverm"
    shell:
        """
        python {WORKDIR}/scripts/merge_coverage.py --files {input} --output {output}
        """