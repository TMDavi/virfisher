import pandas as pd
from Bio import SeqIO
import os
import argparse

def filter_annot(fasta_file, annot_file):
    record_list = [record.id for record in SeqIO.parse(fasta_file, "fasta")]
      # Display first 5 records for debugging
    if not record_list:
        raise ValueError("No records found in the fasta file.")
    
    annot_df = pd.read_csv(annot_file, sep="\t")
    annot_df['contig'] = annot_df['gene'].str.rsplit('_', n=1).str[0]
    filt_annot_df = annot_df[annot_df["contig"].isin(record_list)]
    filt_annot_df = filt_annot_df.drop(columns=['contig'])
    
    return filt_annot_df

def main():
    arg_parser = argparse.ArgumentParser(description="Get final annotation from a CSV file.")
    arg_parser.add_argument("--fasta", type=str, help="Path to the input fasta file.")
    arg_parser.add_argument("--annot", type=str, help="Path to the input annotation file.")
    arg_parser.add_argument("--output", type=str, help="Path to the output CSV file.")
    args = arg_parser.parse_args()

    filtered_df = filter_annot(args.fasta, args.annot)
    filtered_df.to_csv(args.output, index=False, sep="\t")

if __name__ == "__main__":
    main()

