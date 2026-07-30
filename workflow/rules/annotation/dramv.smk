rule table_for_dramv:
    input:
        out("{sample}","intermediate","final_phage_sequences.fasta")
    output:
        viral_seqs = out("{sample}", "intermediate", "virsorter2-step2","final-viral-score.tsv"),
        dramv_table = out("{sample}", "intermediate", "virsorter2-step2", "for-dramv","viral-affi-contigs-for-dramv.tab")
    params:
        outdir = out("{sample}", "intermediate", "virsorter2-step2")
    log:
        stdout = out("{sample}", "intermediate", "virsorter2-step2", "log-stdout.txt"),
        stderr = out("{sample}", "intermediate", "virsorter2-step2", "log-stderr.txt")
    conda:
         "viral-id-sop"
    threads: config["resources"]["threads"]
    shell:
        """
        virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv \
        -i {input} -w vs2-pass2 --include-groups dsDNAphage,ssDNA \
        --min-length 1000 --min-score 0.5 -j {threads} all
        """

rule dramv:
    input:
        viral_seqs = out("{sample}", "intermediate", "virsorter2-step2","final-viral-score.tsv"),
        dramv_table = out("{sample}", "intermediate", "virsorter2-step2", "for-dramv","viral-affi-contigs-for-dramv.tab")
    output:
        out("{sample}","final_results", "annotation","dramv","annotations.tsv")
    params:
        outdir=out("{sample}","final_results", "annotation","dramv")
    conda:
        "viral-id-sop"
    shell:
        """
        DRAM-v.py annotate -i {input.viral_seqs} \
        -v {input.dramv_table} -o {params.outdir} \
        --skip_trnascan --threads {threads} --min_contig_size 1000
        """

rule dramv_distill:
    input:
       out("{sample}","final_results", "annotation","dramv","annotations.tsv")
    output:
        out("{sample}","final_results", "annotation","dramv","dramv-distill","amg_summary.tsv")
    params:
        outdir=out("{sample}","final_results", "annotation","dramv","dramv-distill")
    conda:
        "viral-id-sop"
    shell:
        """
        DRAM-v.py distill -i dramv-annotate/annotations.tsv -o dramv-distill
        """