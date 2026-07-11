import pandas as pd
import argparse
import re
from crass_extract import parse_hmm, parse_domout, extract_genes

def filter_results(input_file, output_file):
    """
    Filters the results based on e-value and score.
    """
    print("Filtering results...")
    df = pd.read_csv(input_file, sep="\t")
    df = df[(df['evalue'] < 1e-10) & (df['score'] > 70) & (df['ali_len'] > 50)]
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
    filtered = gene_counts[gene_counts['gene_hits'] > 5]
    filtered = filtered.sort_values(by='gene_hits', ascending=False)

    if output_filtered_contigs_path=="":
        pass
    else:
        filtered.to_csv(output_filtered_contigs_path, sep="\t", index=False)

    filtered = filtered['contig_id'].tolist()
    print(f"Filtered contigs with more than 5 gene hit: {len(filtered)}")
    return filtered

def main():
    parser = argparse.ArgumentParser(description="Extract putative Crassphage sequences from HMMER results.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file with contigs.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads for HMMER search.")
    parser.add_argument("-o", "--output_dir", required=True, help="Output folder for filtered results.")
    args = parser.parse_args()

    ncldv_db = "/MP_Data/Dados_David/viral_analysis/viral_HMMs/giantDB/NCLDV.hmm"

    parse_hmm(args.input, args.threads, ncldv_db, f"{args.output_dir}/hmmer_results.domout")
    parse_domout(f"{args.output_dir}/hmmer_results.domout", f"{args.output_dir}/parsed_results.txt")
    filter_results(f"{args.output_dir}/parsed_results.txt", f"{args.output_dir}/filtered_results.txt")
    extract_genes(f"{args.output_dir}/filtered_results.txt", 
                  args.input, 
                  count_gene_hits(f"{args.output_dir}/filtered_results.txt"),
                  f"{args.output_dir}/extracted_genes.faa")
    count_gene_hits(f"{args.output_dir}/filtered_results.txt", f"{args.output_dir}/putative_ncldv_contig_list.txt")