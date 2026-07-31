WORKDIR = config["WORKDIR"]
SAMPLES = config["samples"]

rule bwaDB:
    input:
        scaffolds = out("per_sample_results", "{sample}", "intermediate","metaspades", "scaffolds.fasta")
    output:
        amb=out("per_sample_results", "{sample}", "intermediate", "read_alignment", "index", "mappingDB.amb"),
        ann=out("per_sample_results", "{sample}", "intermediate", "read_alignment", "index", "mappingDB.ann"),
        bwt=out("per_sample_results", "{sample}", "intermediate", "read_alignment", "index", "mappingDB.bwt"),
        pac=out("per_sample_results", "{sample}", "intermediate", "read_alignment", "index", "mappingDB.pac"),
        sa=out("per_sample_results", "{sample}", "intermediate", "read_alignment", "index", "mappingDB.sa")
    params:
        index_prefix=out("per_sample_results", "{sample}", "intermediate", "read_alignment", "index","mappingDB")
    shell:
        """
        bwa index -p {params.index_prefix} {input.scaffolds} -a bwtsw
        """

rule bwa_alignment:
    input:
        amb=out("per_sample_results", "{assembly}", "intermediate", "read_alignment", "index", "mappingDB.amb"),
        ann=out("per_sample_results", "{assembly}", "intermediate", "read_alignment", "index", "mappingDB.ann"),
        bwt=out("per_sample_results", "{assembly}", "intermediate", "read_alignment", "index", "mappingDB.bwt"),
        pac=out("per_sample_results", "{assembly}", "intermediate", "read_alignment", "index", "mappingDB.pac"),
        sa=out("per_sample_results", "{assembly}", "intermediate", "read_alignment", "index", "mappingDB.sa"),

        forward=out("per_sample_results", "{reads}", "intermediate", "fastp", "{reads}_R1_trimmed.fastq.gz"),
        reverseR=out("per_sample_results", "{reads}", "intermediate", "fastp", "{reads}_R2_trimmed.fastq.gz"),

    output:
        bam=out("per_sample_results", "{assembly}", "intermediate", "read_alignment","{reads}_vs_{assembly}_sorted.bam")

    log:
        stderr=out("per_sample_results", "{assembly}", "intermediate", "read_alignment","logs", "{reads}_vs_{assembly}_bwa.stderr.log")

    benchmark:
        out("per_sample_results", "{assembly}", "intermediate", "read_alignment","benchmarks", "{reads}_vs_{assembly}_bwa.txt")

    threads:
        config["resources"]["threads"]

    params:
        index_prefix=out("per_sample_results", "{assembly}", "intermediate", "read_alignment",
                         "index", "mappingDB")

    conda:
        "coverm"

    shell:
        """
        bwa mem -t {threads} {params.index_prefix} {input.forward} {input.reverseR} 2> {log.stderr} | \
            samtools sort -@ {threads} -o {output.bam}
        """

rule metabat_coverage:
    input:
        bams=lambda wc: expand(
            out("per_sample_results", wc.assembly, "intermediate","read_alignment","{read}_vs_{assembly}_sorted.bam"),
            read=config["samples"],
            assembly=wc.assembly
        )
    output:
        out("per_sample_results", "{assembly}", "final_results", "read_alignment","{assembly}_metabat.coverage")

    conda:
        "coverm"

    shell:
        """
        coverm contig -b {input.bams} --methods metabat -o {output}
        """

rule scaffold_coverage:
    input:
        bams=lambda wc: expand(
            out("per_sample_results", wc.assembly,"intermediate","read_alignment","{read}_vs_{assembly}_sorted.bam"),
            read=config["samples"],
            assembly=wc.assembly
        )

    output:
        out("per_sample_results", "{assembly}", "final_results", "read_alignment","{assembly}.coverage")

    conda:
        "coverm"

    shell:
        """
        coverm contig -b {input.bams} --methods count covered_fraction -o {output}
        """