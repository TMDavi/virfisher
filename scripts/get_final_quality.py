import pandas as pd
from Bio import SeqIO
import os
import argparse

def filter_quality(fasta_file, qual_file):
    record_list = [record.id for record in SeqIO.parse(fasta_file, "fasta")]
      # Display first 5 records for debugging
    if not record_list:
        raise ValueError("No records found in the fasta file.")
    
    quality_df = pd.read_csv(qual_file, sep="\t")
    quality_df['contig_id'] = quality_df['contig_id'].str.split(r"\|\|").str[0]
    filt_quality_df = quality_df[quality_df["contig_id"].isin(record_list)]
    
    return filt_quality_df

def main():
    arg_parser = argparse.ArgumentParser(description="Get final quality from a CSV file.")
    arg_parser.add_argument("--fasta", type=str, help="Path to the input fasta file.")
    arg_parser.add_argument("--qual", type=str, help="Path to the input quality file.")
    arg_parser.add_argument("--output", type=str, help="Path to the output CSV file.")
    args = arg_parser.parse_args()

    filtered_df = filter_quality(args.fasta, args.qual)
    filtered_df.to_csv(args.output, index=False, sep="\t")

if __name__ == "__main__":
    main()