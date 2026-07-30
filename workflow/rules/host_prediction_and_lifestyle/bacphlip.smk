rule bacphlip:
    input:
        out("{sample}","intermediate", "final_phage_sequences.fasta")
    output:
        directory(out("{sample}","final_results","life_style","bacphlip"))
    conda:
        "bacphlip"
    params:
        outdir=out("{sample}","intermediate", "final_phage_sequences.fasta.BACPHLIP_DIR"),
        bacfile=out("{sample}","intermediate", "final_phage_sequences.fasta.bacphlip"),
        hmmfile=out("{sample}","intermediate", "final_phage_sequences.fasta.hmmsearch.tsv")
    shell:
        """
        bacphlip -i {input} --multi_fasta -f

        rm -r {params.outdir}

        mkdir {output}

        mv {params.bacfile} {output}
        mv {params.hmmfile} {output}
        
        """