
rule rename_contigs:
    input:  
        assembly=out("{sample}","intermediate", "viral_predicted_scaffolds_second_step.fasta")
    output:
        temp(out("{sample}", "intermediate", "dereplication", "{sample}_scaffolds.fasta"))
    params:
        sample_name=lambda wc: wc.sample
    shell:
        """
        python {WORKDIR}/scripts/rename_contigs.py -i {input.assembly} -o {output} -p {params.sample_name}
        """
rule merge_contigs:
    input:
        expand(
            out("{sample}", "intermediate", "dereplication", "{sample}_scaffolds.fasta"),
            sample=list(config["samples"].keys())
        )
    output:
        temp(out("dereplicated", "merged_scaffolds.fasta"))
    shell:
        """
        cat {input} > {output}
        """

rule makeblastdb:
    input:
        out("dereplicated", "merged_scaffolds.fasta")
    output:
        ndb = out("dereplicated","blastdb","all_viral_contigs.ndb"),
        nhr = out("dereplicated","blastdb","all_viral_contigs.nhr"),
        nin = out("dereplicated","blastdb","all_viral_contigs.nin"),
        no = out("dereplicated","blastdb","all_viral_contigs.not"),
        nsq = out("dereplicated","blastdb","all_viral_contigs.nsq"),
        ntf = out("dereplicated","blastdb","all_viral_contigs.ntf"),
        nto = out("dereplicated","blastdb","all_viral_contigs.nto")
    params:
        db_name = out("dereplicated","blastdb","all_viral_contigs")
    conda:
        "viral-id-sop"
    shell:
        """
        makeblastdb -in {input} -dbtype nucl -out {params.db_name}
        """

rule blastn:
    input:
        query = out("dereplicated", "merged_scaffolds.fasta"),
        ndb = out("dereplicated","blastdb","all_viral_contigs.ndb"),
        nhr = out("dereplicated","blastdb","all_viral_contigs.nhr"),
        nin = out("dereplicated","blastdb","all_viral_contigs.nin"),
        no = out("dereplicated","blastdb","all_viral_contigs.not"),
        nsq = out("dereplicated","blastdb","all_viral_contigs.nsq"),
        ntf = out("dereplicated","blastdb","all_viral_contigs.ntf"),
        nto = out("dereplicated","blastdb","all_viral_contigs.nto")

    output:
        out("dereplicated","all_viral_contigs_blastn.tsv")

    params:
        db = out("dereplicated","blastdb","all_viral_contigs")
    conda:
        "viral-id-sop"
    threads:config["resources"]["threads"]
    shell:
        """
        blastn -query {input.query} -db {params.db} -outfmt '6 std qlen slen' -max_target_seqs 10000 -out {output} -num_threads {threads}
        """

rule anicalc:
    input:
        out("dereplicated","all_viral_contigs_blastn.tsv")
    output:
        out("dereplicated","all_viral_contigs_ani.tsv")
 
    shell:
        """
        python programs/clustering_scripts/anicalc.py -i {input} -o {output}
        """

rule aniclust:
    input:
        fasta = out("dereplicated", "merged_scaffolds.fasta"),
        ani = out("dereplicated","all_viral_contigs_ani.tsv")
    output:
        out("dereplicated","all_viral_contigs_clusters.tsv")
    shell:
        """
        python programs/clustering_scripts/aniclust.py --fna {input.fasta} --ani {input.ani} --out {output} --min_ani 95 --min_tcov 85 --min_qcov 0
        """

rule filter_miuvigs:
    input:
        fasta= out("dereplicated", "merged_scaffolds.fasta"),
        derep_miuvigs_list= out("dereplicated","all_viral_contigs_clusters.tsv")
    output:
        fasta = out("dereplicated","derep_miuvigs.fasta"),
        ids = out("dereplicated","miuvig_ids.txt")
    shell:
        """
        cut -f 1 {input.derep_miuvigs_list} > {output.ids}

        seqkit grep -f {output.ids} {input.fasta} > {output.fasta}
        """

rule redistribute_miuvigs:
    input:
        out("dereplicated","derep_miuvigs.fasta")
    output:
        expand(out("{sample}","final_results","dereplicated","miuvigs.fasta"), sample=config['samples'].keys())
    params:
        samplelist = list(config["samples"].keys()),
        outdir=out("")
    shell:
        """
        python {WORKDIR}/scripts/distribute_contigs.py -m {input} -s {params.samplelist} -o {params.outdir}
        """