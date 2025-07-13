import pandas as pd
import re
import os
import argparse
from crass_extract import parse_hmm, parse_domout, extract_genes 


def filter_virophage_results(input_file, threshold_file, output_file):
    """
    Filters the results based on e-value and score.
    """
    thresholds = pd.read_csv(threshold_file, sep="\t")
    print("Filtering results...")
    df = pd.read_csv(input_file, sep="\t")
    df = df[(df['evalue'] < 1e-10) & (df['score'] > 70)]

    mask = []  
    for index, row in df.iterrows():
        gene = row['qname'].split('_')[0]
        if gene in thresholds['Gene'].values:
            cutoff = thresholds.loc[thresholds['Gene'] == gene, 'ali_cut'].values[0]

            mask.append(row['ali_len'] >= cutoff)
        else:
            mask.append(True)
    df = df[mask]  

    df['contig_name'] = df['contig_id'].astype(str).str.replace(r'_\d+$', '', regex=True)
    df = df.sort_values(by='score', ascending=False)
    df = df.drop_duplicates(subset=['contig_name', 'qname'], keep='first')
    df = df.drop('contig_name',axis=1)
    df.to_csv(output_file, sep="\t", index=False)

def count_virophage_gene_hits(input_file, output_filtered_contigs_path=""):
    """
    Counts the number of gene hits per contig and filters those with more than 4 hits.
    """
    print("Counting gene hits...")
    df = pd.read_csv(input_file, sep="\t")
    df['contig_id'] = df['contig_id'].astype(str).apply(lambda x: re.sub(r'_\d+$', '', x))
    
    genes = ['MCP', 'Penton', 'ATPase', 'PRO']
    for gene in genes:
        df[gene] = df['qname'].astype(str).str.contains(gene, case=False).astype(int)

    gene_counts_per_contig = df.groupby('contig_id')[genes].sum().reset_index()
    gene_counts_per_contig = gene_counts_per_contig[gene_counts_per_contig['MCP'] > 0]

    if output_filtered_contigs_path:
        gene_counts_per_contig.to_csv(output_filtered_contigs_path, sep="\t", index=False)

    filtered = gene_counts_per_contig['contig_id'].tolist()
    print(f"Filtered contigs with at least MCP hit: {len(filtered)}")
    return filtered

def filter_plv_results(input_file, output_file):
    """
    Filters the results based on e-value and score for PLV.
    """

    virophage_contigs = count_virophage_gene_hits(input_file)

    print("Filtering PLV results...")
    df = pd.read_csv(input_file, sep="\t")
    df = df[(df['evalue'] < 1e-10) & (df['score'] > 70) & (df['ali_len'] > 50)]
    df['contig_name'] = df['contig_id'].astype(str).str.replace(r'_\d+$', '', regex=True)
    df = df[~df['contig_name'].isin(virophage_contigs)]#Removing any conitgs that have been previously classified as virophages
    df = df.sort_values(by='score', ascending=False)
    df = df.drop_duplicates(subset=['contig_name', 'qname'], keep='first')
    df = df.drop('contig_name', axis=1)
    df.to_csv(output_file, sep="\t", index=False)

def count_plv_gene_hits(input_file, output_filtered_contigs_path=""):
    print("Counting PLV gene hits...") 
    df = pd.read_csv(input_file, sep="\t")
    df['contig_id'] = df['contig_id'].astype(str).apply(lambda x: re.sub(r'_\d+$', '', x))
    gene_counts = df.groupby('contig_id').size().reset_index(name='gene_hits')
    filtered = gene_counts[gene_counts['gene_hits'] > 0]
    filtered = filtered.sort_values(by='gene_hits', ascending=False)

    if output_filtered_contigs_path:
        filtered.to_csv(output_filtered_contigs_path, sep="\t", index=False) 

    filtered_list = filtered['contig_id'].tolist()
    print(f"Filtered contigs with more than 1 gene hit: {len(filtered)}")
    return filtered_list


def main():
    parser = argparse.ArgumentParser(description="Extract virophage sequences from a FASTA file based on HMMER results.")
    parser.add_argument("-i","--input", help="Input FASTA file containing protein sequences.")
    parser.add_argument("-t","--threads", help="Number of threads.")
    parser.add_argument("-o","--output_dir", help="Output directory to save filtered virophage sequences.")
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    threshold_file = "/MP_Data/Dados_David/viral_analysis/viral_HMMs/virophageDB/cutoffs.tsv"
    if not os.path.exists(threshold_file):
        raise FileNotFoundError(f"Threshold file {threshold_file} does not exist.")

    #For virophage identification

    virophage_db = "/MP_Data/Dados_David/viral_analysis/viral_HMMs/virophageDB/virophages_All_markers.hmm"
    if not os.path.exists(virophage_db):
        raise FileNotFoundError(f"Virophage HMM database {virophage_db} does not exist.")

    parse_hmm(args.input, args.threads, virophage_db, f"{args.output_dir}/virophage_hmmer_results.domout")
    parse_domout(f"{args.output_dir}/virophage_hmmer_results.domout", f"{args.output_dir}/virophage_parsed_results.txt")
    filter_virophage_results(f"{args.output_dir}/virophage_parsed_results.txt", threshold_file, f"{args.output_dir}/virophage_filtered_results.txt")
    extract_genes(f"{args.output_dir}/virophage_filtered_results.txt", 
                  args.input, 
                  count_virophage_gene_hits(f"{args.output_dir}/virophage_filtered_results.txt"),
                  f"{args.output_dir}/virophage_extracted_genes.faa")
    count_virophage_gene_hits(f"{args.output_dir}/virophage_filtered_results.txt", f"{args.output_dir}/virophage_filtered_contigs.txt")

    #For plv identification
    plv_db = "/MP_Data/Dados_David/viral_analysis/viral_HMMs/virophageDB/PLV.hmm"
    if not os.path.exists(plv_db):
        raise FileNotFoundError(f"PLV HMM database {plv_db} does not exist.")
    parse_hmm(args.input, args.threads, plv_db, f"{args.output_dir}/plv_hmmer_results.domout")
    parse_domout(f"{args.output_dir}/plv_hmmer_results.domout", f"{args.output_dir}/plv_parsed_results.txt")
    filter_plv_results(f"{args.output_dir}/plv_parsed_results.txt", f"{args.output_dir}/plv_filtered_results.txt")
    extract_genes(f"{args.output_dir}/plv_filtered_results.txt", 
                  args.input,
                  count_plv_gene_hits(f"{args.output_dir}/plv_filtered_results.txt"), 
                  f"{args.output_dir}/plv_extracted_genes.faa")
    count_plv_gene_hits(f"{args.output_dir}/plv_filtered_results.txt", f"{args.output_dir}/plv_filtered_contigs.txt")

if __name__ == "__main__":
    main()