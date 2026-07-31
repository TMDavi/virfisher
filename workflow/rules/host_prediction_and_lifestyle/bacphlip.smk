rule bacphlip:
    input:
        out("per_sample_results", "{sample}","intermediate", "final_phage_sequences.fasta")
    output:
        directory(out("per_sample_results", "{sample}","final_results","life_style","bacphlip"))
    conda:
        "bacphlip"
    params:
        outdir=out("per_sample_results", "{sample}","intermediate", "final_phage_sequences.fasta.BACPHLIP_DIR"),
        bacfile=out("per_sample_results", "{sample}","intermediate", "final_phage_sequences.fasta.bacphlip"),
        hmmfile=out("per_sample_results", "{sample}","intermediate", "final_phage_sequences.fasta.hmmsearch.tsv")
    shell:
        """
        bacphlip -i {input} --multi_fasta -f

        rm -r {params.outdir}

        mkdir {output}

        mv {params.bacfile} {output}
        mv {params.hmmfile} {output}
        
        """