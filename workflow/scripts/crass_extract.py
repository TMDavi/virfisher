import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import re
import argparse


def parse_hmm(input_file, threads, dbname, hmmerfile):
    """
    Parses the HMM file and returns a list of sequences.
    """
    print("Performing HMM search...")
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file {input_file} does not exist.")
    if not os.path.exists(dbname):
        raise FileNotFoundError(f"CrassDB file {dbname} does not exist.")
    cmd = cmd = f"hmmsearch --cpu {threads} -E 1e-3 --domtblout {hmmerfile} {dbname} {input_file}"
    os.system(cmd)

def parse_domout(hmmerfile, output_file):
    """
    Parses the HMMER domout file and writes the results to a new file.
    """
    print("Parsing HMMER domout file...")
    with open(hmmerfile, "r") as f, open(output_file,"w") as parsed_file:

        parsed_file.write("contig_id\ttlen\tqname\tali_from\tali_to\tali_len\tqlen\tscore\tevalue\n")

        domout = f.readlines()
        lines = [line.strip().split() for line in domout if line.strip()]
        for line in lines:
            if line[0].startswith("#"):
                continue
            contig_id = line[0]
            tlen = line[2]
            qname = line[3]
            qlen= line[5]
            evalue=line[6]
            score = float(line[7])
            ali_from = int(line[17])
            ali_to = int(line[18])
            ali_len = abs(ali_to - ali_from) + 1
            parsed_file.write(f"{contig_id}\t{tlen}\t{qname}\t{ali_from}\t{ali_to}\t{ali_len}\t{qlen}\t{score}\t{evalue}\n")

def filter_results(input_file, output_file):
    """
    Filters the results based on e-value and score.
    """
    print("Filtering results...")
    df = pd.read_csv(input_file, sep="\t")
    df = df[(df['evalue'] < 1e-10) & (df['score'] > 70) & (df['ali_len'] > 30)]
    df['contig_name'] = df['contig_id'].astype(str).str.replace(r'_\d+$', '', regex=True)
    df = df.sort_values(by='score', ascending=False)
    df = df.drop_duplicates(subset=['contig_name', 'qname'], keep='first')
    df = df.drop('contig_name',axis=1)
    df.to_csv(output_file, sep="\t", index=False)

def count_gene_hits(input_file, output_filtered_contigs_path=""):
    """
    Counts the number of gene hits per contig and filters those with more than 4 hits.
    """
    print("Counting gene hits...")
    df = pd.read_csv(input_file, sep="\t")
    df['contig_id'] = df['contig_id'].astype(str).apply(lambda x: re.sub(r'_\d+$', '', x))
    gene_counts = df.groupby('contig_id').size().reset_index(name='gene_hits')
    filtered = gene_counts[gene_counts['gene_hits'] > 1]
    filtered = filtered.sort_values(by='gene_hits', ascending=False)

    if output_filtered_contigs_path=="":
        pass
    else:
        filtered.to_csv(output_filtered_contigs_path, sep="\t", index=False)

    filtered = filtered['contig_id'].tolist()
    print(f"Filtered contigs with more than 1 gene hit: {len(filtered)}")
    return filtered


def extract_genes(input_file, fasta_file, final_contig_list, output_fasta):
    seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    extracted_records = []

    #final_contig_list = count_gene_hits(input_file)
    df = pd.read_csv(input_file, sep="\t")
    df['contig_name'] = df['contig_id'].astype(str).str.replace(r'_\d+$', '', regex=True)
    df = df[df['contig_name'].isin(final_contig_list)]
    df = df.drop('contig_name', axis=1)

    for index, row in df.iterrows():
        contig_id = row['contig_id']
        query_name = row['qname']
        ali_from = int(row['ali_from'])
        ali_to = int(row['ali_to'])
        score = float(row['score'])
        evalue = float(row['evalue'])

        if contig_id not in seq_dict:
            print(f"Warning: {contig_id} not found in FASTA.")
            continue

        # Handle reverse strands if needed
        start = min(ali_from, ali_to) - 1  
        end = max(ali_from, ali_to)

        sequence = seq_dict[contig_id].seq[start:end]
        if ali_from > ali_to:
            sequence = sequence.reverse_complement()

        new_id = f"{contig_id}_{query_name}_{start+1}_{end}"
        record = SeqRecord(sequence, id=new_id, description=f"E-value={evalue}; Score={score}")
        extracted_records.append(record)

    SeqIO.write(extracted_records, output_fasta, "fasta")



def main():
    parser = argparse.ArgumentParser(description="Extract putative Crassphage sequences from HMMER results.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file with contigs.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads for HMMER search.")
    parser.add_argument("-o", "--output_dir", required=True, help="Output folder for filtered results.")
    args = parser.parse_args()

    crassdb = "../viral_HMMs/crassDB/crass_hallmark.hmm"

    parse_hmm(args.input, args.threads, crassdb, f"{args.output_dir}/hmmer_results.domout")
    parse_domout(f"{args.output_dir}/hmmer_results.domout", f"{args.output_dir}/parsed_results.txt")
    filter_results(f"{args.output_dir}/parsed_results.txt", f"{args.output_dir}/filtered_results.txt")
    extract_genes(f"{args.output_dir}/filtered_results.txt", 
                  args.input, 
                  count_gene_hits(f"{args.output_dir}/filtered_results.txt"),
                  f"{args.output_dir}/extracted_genes.faa")
    count_gene_hits(f"{args.output_dir}/filtered_results.txt", f"{args.output_dir}/putative_crass_contig_list.txt")

if __name__ == "__main__":
    main()
    


  