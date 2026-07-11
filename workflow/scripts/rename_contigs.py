from Bio import SeqIO
import os
import argparse

def rename_contigs(input_file, output_file, sample_name):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            record.id = f"{sample_name}_{record.id}"
            SeqIO.write(record, outfile, "fasta")
    return output_file

def main():
    parser = argparse.ArgumentParser(description="Rename contigs in a FASTA file with a sample name prefix.")
    parser.add_argument("-i","--input_file", type=str, help="Input FASTA file containing contigs.")
    parser.add_argument("-o","--output_file", type=str, help="Output FASTA file with renamed contigs.")
    parser.add_argument("-p","--prefix", type=str, help="Sample name to prefix to each contig ID.")
    
    args = parser.parse_args()
    
    rename_contigs(args.input_file, args.output_file, args.prefix)

if main() == "__main__":
    main()