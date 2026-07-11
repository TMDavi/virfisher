from Bio import SeqIO
import os
import pandas as pd
import argparse

def treat_checkvID(viruses_file, proviruses_file, combined_file):
    seen_ids = set()
    unique_records = []

    # Parse virus sequences
    for record in SeqIO.parse(viruses_file, "fasta"):
        identifier = record.id.split(r"\|\|")[0]
        if identifier in seen_ids:
            continue  # Skip duplicates
        seen_ids.add(identifier)
        record.id = identifier
        unique_records.append(record)

    # Parse provirus sequences
    for record in SeqIO.parse(proviruses_file, "fasta"):
        if record.id.endswith("_1"):
            identifier = "_".join(record.id.split("_")[:-1])
        else:
            identifier = record.id
        identifier = identifier.split(r"\|\|")[0]
        if identifier in seen_ids:
            continue  # Skip duplicates
        seen_ids.add(identifier)
        record.id = identifier
        unique_records.append(record)

    # Write only unique records
    SeqIO.write(unique_records, combined_file, "fasta")
    return combined_file

def main():
    parser = argparse.ArgumentParser(description="Combine CheckV IDs from viruses and proviruses FASTA files.")
    parser.add_argument("--virus", help="Path to the viruses FASTA file.")
    parser.add_argument("--provirus", help="Path to the proviruses FASTA file.")
    parser.add_argument("--output", help="Path to the output combined FASTA file.")

    args = parser.parse_args()

    combined_file = treat_checkvID(args.virus, args.provirus, args.output)
    print(f"✅ Combined FASTA file created at: {combined_file}")

if __name__ == "__main__":
    main()
    # Example usage:
    # python combine_checkv.py --virus viruses.fasta --provirus proviruses.fasta --output combined.fasta