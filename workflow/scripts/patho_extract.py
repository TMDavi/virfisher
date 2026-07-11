import pandas as pd
import os
import argparse
import re
from crass_extract import parse_hmm, parse_domout, extract_genes

def get_set_repeated_values(list_of_sets):
    """
    Gets values that repeat inside a list of sets
    """

    if not list_of_sets:
        return set()

    repeated_values = set()

    for i in range(len(list_of_sets)):
        for j in range(i + 1, len(list_of_sets)):
            repeated_values.update(list_of_sets[i].intersection(list_of_sets[j]))
    
    return repeated_values

    

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


def separate_contigs_by_taxonomy(input_file, taxonomy_file):
    
    #Separates contigs by taxonomy and saves them in separate lists.
   
    adeno_list = []
    herpes_list = []
    papiloma_list = []
    polyoma_list = []

    taxonomy_df = pd.read_csv(taxonomy_file, sep="\t")
    gene_to_taxonomy = dict(zip(taxonomy_df['Gene'], taxonomy_df['Taxonomy']))

    print("Separating contigs by taxonomy...")
    df = pd.read_csv(input_file, sep="\t")
    df['contig_id'] = df['contig_id'].astype(str).apply(lambda x: re.sub(r'_\d+$', '', x))

    for index, row in df.iterrows():
        gene = row['qname']
        contig_id = row['contig_id']
        taxon = gene_to_taxonomy(gene)
        if taxon == 'Adenoviridae':
            adeno_list.append(contig_id)
        elif taxon == 'Herpesviridae':
            herpes_list.append(contig_id)
        elif taxon == 'Papillomaviridae':
            papiloma_list.append(contig_id)
        elif taxon == 'Polyomaviridae':
            polyoma_list.append(contig_id)

    adeno_set = set(adeno_list)
    herpes_set = set(herpes_list)
    papiloma_set = set(papiloma_list)
    polyoma_set = set(polyoma_list)

    repeated = get_set_repeated_values([adeno_set, herpes_set, papiloma_set, polyoma_set])


def main():
    parser = argparse.ArgumentParser(description="Extract virophage sequences from a FASTA file based on HMMER results.")
    parser.add_argument("-i","--input", help="Input FASTA file containing protein sequences.")
    parser.add_argument("-t","--threads", help="Number of threads.")
    parser.add_argument("-o","--output_dir", help="Output directory to save filtered virophage sequences.")
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    tax_file = "viral_HMMs/pathoDB/tax_assignment.tsv"
    patho_db = "viral_HMMs/pathoDB/patho_hmm_db.hmm"

# Example Usage:
set1 = {1, 2, 3}
set2 = {3, 4, 5, 6}
set3 = {6, 7, 8, 1}
set4 = {9, 10}

all_my_sets = [set1, set2, set3, set4]
result = get_set_repeated_values(all_my_sets)
print(f"repeated values: {result}")
print(f"set 1 unique values {set1 - result}")