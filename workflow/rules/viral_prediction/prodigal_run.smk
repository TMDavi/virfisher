rule separate_contigs:
    input:
        assembly=out("{sample}", "intermediate", "metaspades", "scaffolds.fasta")
    output:
        contigs10kb=out("{sample}", "intermediate", "NCLDV_id", "input.min10kb.fasta"),
        contigdir=directory(out("{sample}", "intermediate", "NCLDV_id", "contigs"))
    params:
        outdir=out("{sample}", "intermediate", "NCLDV_id", "contigs")
    shell:
        r"""
        mkdir -p "{params.outdir}"

        seqkit seq -m 10000 "{input.assembly}" > "{output.contigs10kb}"

        seqkit faidx "{output.contigs10kb}"

        seqkit seq -n "{output.contigs10kb}" | while read id; do
            short_id=$(echo "$id" | sed -E 's/^(NODE_[0-9]+).*/\1/')
            seqkit faidx "{output.contigs10kb}" "$id" > "{params.outdir}/${{short_id}}.fasta"
        done
        """

rule prodigal_gv:
    input:
        contigdir = directory(out("{sample}", "intermediate", "NCLDV_id", "contigs"))
    output:
        annotdir = directory(out("{sample}", "intermediate", "NCLDV_id", "prodigal"))
    conda:
        "ncdlv_msearch"
    threads: config["resources"]["threads"]
    shell:
        r"""
        mkdir -p "{output.annotdir}"

        for file in "{input.contigdir}"/*.fasta; do
            base=$(basename "$file" .fasta)

            python {WORKDIR}/scripts/parallel-prodigal-gv.py -t {threads} -q -i "$file" -a "{output.annotdir}/${{base}}.faa"
        done
        """

rule ncldv_msearch:
    input:
        out("{sample}", "intermediate", "NCLDV_id", "prodigal")
    output:
        full=out("{sample}", "intermediate", "NCLDV_id", "results","{sample}.full_output.txt")
    conda:
        "ncdlv_msearch"
    params:
        samplename="{sample}",
        outdir=out("{sample}", "intermediate", "NCLDV_id", "results")
    threads:
        config["resources"]["threads"]
    shell:
        """
        python programs/ncldv_markersearch/ncldv_markersearch.py -i {input} -a \
        -m GVOGm0003,GVOGm0890,GVOGm0760,GVOGm0461,GVOGm0172,GVOGm0054,GVOGm0023,GVOGm0013,GVOGm0022 \
        -t {threads} -n {params.samplename}

        cp {params.samplename}* {params.outdir}
        rm {params.samplename}*
        """