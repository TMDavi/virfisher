from Bio import SeqIO
import os
from collections import defaultdict
import argparse

def extract_miuvig(miuvig_file, samplelist, outdir):
    sample_set = set(samplelist)
    sample_records = defaultdict(list)

    for record in SeqIO.parse(miuvig_file, "fasta"):
        samplename = record.id.split("_NODE", 1)[0]
        if samplename in sample_set:
            sample_records[samplename].append(record)
    
    for sample, records in sample_records.items():
        if not os.path.exists(f"{outdir}/{sample}/final_results/dereplicated"):
            os.makedirs(f"{outdir}/{sample}/final_results/dereplicated")
        outfile = f"{outdir}/{sample}/final_results/dereplicated/miuvigs.fasta"
        SeqIO.write(records, outfile, "fasta")

def main():
    parser = argparse.ArgumentParser("Redristibute the clustered MiUviGs")
    parser.add_argument("-m", "--miuvigs", help="The clustered miuvig file")
    parser.add_argument("-s", "--samples",nargs="+", help="The list of samples to redistribute the miuvigs")
    parser.add_argument("-o","--outdir", help="The output directory")
    
    args = parser.parse_args()

    extract_miuvig(args.miuvigs, args.samples, args.outdir)

if __name__ == "__main__":
    main()
