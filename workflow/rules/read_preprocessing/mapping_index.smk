WORKDIR = config["WORKDIR"]
SAMPLES = config["samples"]

rule bwaDB:
    input:
        scaffolds = out("{sample}", "intermediate","metaspades", "scaffolds.fasta")
    output:
        amb=out("{sample}", "intermediate", "read_alignment", "index", "mappingDB.amb"),
        ann=out("{sample}", "intermediate", "read_alignment", "index", "mappingDB.ann"),
        bwt=out("{sample}", "intermediate", "read_alignment", "index", "mappingDB.bwt"),
        pac=out("{sample}", "intermediate", "read_alignment", "index", "mappingDB.pac"),
        sa=out("{sample}", "intermediate", "read_alignment", "index", "mappingDB.sa")
    params:
        index_prefix=out("{sample}", "intermediate", "read_alignment", "index","mappingDB")
    shell:
        """
        bwa index -p {params.index_prefix} {input.scaffolds} -a bwtsw
        """

rule bwa_alignment:
    input:
        amb=out("{assembly}", "intermediate", "read_alignment", "index", "mappingDB.amb"),
        ann=out("{assembly}", "intermediate", "read_alignment", "index", "mappingDB.ann"),
        bwt=out("{assembly}", "intermediate", "read_alignment", "index", "mappingDB.bwt"),
        pac=out("{assembly}", "intermediate", "read_alignment", "index", "mappingDB.pac"),
        sa=out("{assembly}", "intermediate", "read_alignment", "index", "mappingDB.sa"),

        forward=out("{reads}", "intermediate", "fastp", "{reads}_R1_trimmed.fastq.gz"),
        reverseR=out("{reads}", "intermediate", "fastp", "{reads}_R2_trimmed.fastq.gz"),

    output:
        bam=out("{assembly}", "intermediate", "read_alignment","{reads}_vs_{assembly}_sorted.bam")

    log:
        stderr=out("{assembly}", "intermediate", "read_alignment","logs", "{reads}_vs_{assembly}_bwa.stderr.log")

    benchmark:
        out("{assembly}", "intermediate", "read_alignment","benchmarks", "{reads}_vs_{assembly}_bwa.txt")

    threads:
        config["resources"]["threads"]

    params:
        index_prefix=out("{assembly}", "intermediate", "read_alignment",
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
            out(wc.assembly, "intermediate","read_alignment","{read}_vs_{assembly}_sorted.bam"),
            read=config["samples"],
            assembly=wc.assembly
        )
    output:
        out("{assembly}", "final_results", "read_alignment","{assembly}_metabat.coverage")

    conda:
        "coverm"

    shell:
        """
        coverm contig -b {input.bams} --methods metabat -o {output}
        """

rule scaffold_coverage:
    input:
        bams=lambda wc: expand(
            out(wc.assembly,"intermediate","read_alignment","{read}_vs_{assembly}_sorted.bam"),
            read=config["samples"],
            assembly=wc.assembly
        )

    output:
        out("{assembly}", "final_results", "read_alignment","{assembly}.coverage")

    conda:
        "coverm"

    shell:
        """
        coverm contig -b {input.bams} --methods count covered_fraction -o {output}
        """