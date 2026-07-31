rule vibrant:
    input:
        out("per_sample_results", "{sample}","intermediate", "final_phage_sequences.fasta")
    output:
        lyso = out("per_sample_results", "{sample}","final_results","life_style","vibrant_results", "final_phage_sequences.phages_lysogenic.fna"),
        lytic = out("per_sample_results", "{sample}","final_results","life_style","vibrant_results", "final_phage_sequences.phages_lytic.fna")
    conda:
        "vibrant"
    params:
        outdir=out("per_sample_results", "{sample}","intermediate","vibrant"),
        outdirtoremove=out("per_sample_results", "{sample}","intermediate","vibrant","VIBRANT_final_phage_sequences"),
        lyso = out("per_sample_results", "{sample}","intermediate","vibrant","VIBRANT_final_phage_sequences","VIBRANT_phages_final_phage_sequences","final_phage_sequences.phages_lysogenic.fna"),
        lytic = out("per_sample_results", "{sample}","intermediate","vibrant","VIBRANT_final_phage_sequences","VIBRANT_phages_final_phage_sequences","final_phage_sequences.phages_lytic.fna"),
        finaldir=out("per_sample_results", "{sample}","final_results","life_style")
    threads:
        config["resources"]["threads"]
    shell:
        """
        VIBRANT_run.py -i {input} -f nucl -t {threads} -folder {params.outdir}

        mkdir -p {params.finaldir}
        mkdir -p {params.finaldir}/vibrant_results
        cp {params.lyso} {params.finaldir}/vibrant_results
        cp {params.lytic} {params.finaldir}/vibrant_results

        rm -r {params.outdirtoremove}
        """