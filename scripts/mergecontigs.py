from Bio import SeqIO
import argparse
import os

def extract_contigs(fasta_file, file_list): 

    seen_ids = set()

    with open(fasta_file, "w") as mergedcontigs:
        for file in file_list:
            if not os.path.exists(file):
                print(f"File {file} does not exist.")
                continue

            for record in SeqIO.parse(file, "fasta"):
                ide = record.id.split("||")[0]
                if ide not in seen_ids:
                    SeqIO.write(record, mergedcontigs, "fasta")
                    seen_ids.add(ide)
                else:
                    print(f"Duplicate ID skipped: {ide}")

    print(f"\nTotal unique headers written: {len(seen_ids)}")

def main():
    parser = argparse.ArgumentParser(description="combine different fasta files with repeated ids")
    parser.add_argument("--output", help="Output file with unique headers.")
    parser.add_argument("--vs2_file", help="Virsorter2 file")
    parser.add_argument("--dvf_file", help="Deepvirfinder file")

    args = parser.parse_args()

    file_list = [args.vs2_file, args.dvf_file]

    extract_contigs(args.output, file_list)

if __name__ == "__main__":
    main()