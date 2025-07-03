import pandas as pd
import re
import os
import argparse
from Bio import SeqIO

suspicious = ['carbohydrate kinase','carbohydrate-kinase','glycosyltransferase',
                'glycosyl transferase','glycosyl transferaseendonuclease','nucleotide sugar epimerase',
                'nucleotide sugar-epimerase','nucleotide-sugar epimerase','nucleotide-sugar-epimerase',
                'nucleotidyltransferase','nucleotidyl transferase','nucleotidyl-transferase', 
                'plasmid stability','endonuclease','ABC transporter','CRISPR Cas','Sporulation',
                'Two-component system','Secretion system']    

def changeName(s):
    return re.sub(r'_\d+$', '', s) #Remove the trailing _1, _2, etc. from the gene name

def findSuspicious(gene_name):
    """
    Check if the gene name contains any of the suspicious substrings.
    """
    if isinstance(gene_name, str):  # Check if the gene_name is a string
        for substring in suspicious:
            if substring in gene_name:
                return True
    return False    


def getSuspiciosContigs(annotation_file):
    """
    Extracts contig IDs from the annotation file based on suspicious gene names.
    If more than half of the genes in a contig are suspicious, the contig is removed
    """
    df = pd.read_csv(annotation_file, sep='\t')
    
    df['is_suspicious'] = df['annotation_description'].apply(findSuspicious)
    df['gene'] = df['gene'].apply(changeName)
    df = df[df['annotation_description'].notna()]
    df = df[df['plasmid_hallmark'] == 0] #Filter contigs with plasmid hallmark

    #Set df1
    df1 = df['gene']
    df1 = df1.value_counts().rename_axis('gene').reset_index(name='gene_counts')

    df2 = df[['gene','is_suspicious']]
    df2 = df2[df2['is_suspicious'] == True].dropna()
    df2 = df2.groupby(['gene', 'is_suspicious']).size().reset_index(name='sus_count')

    result = pd.merge(df1, df2, how='outer',on='gene')
    result= result[['gene','gene_counts', 'sus_count']]

    to_remove= result[(result['sus_count'] * 2) >= result['gene_counts']]
    #to_remove = to_remove['contig_id']
    #print(to_remove)
    final = result[~result['gene'].isin(to_remove['gene'])]
    #final.to_csv('filtered_genes.csv', index=False)
    #final = final[final['gene_counts'] >= 3]
    final = final['gene'].tolist()
    return final

def extractContigs(fasta_file, annotation_file, output_file):
    
    headers_to_keep = getSuspiciosContigs(annotation_file)
    headers_to_keep = set(headers_to_keep)  # Convert to set for faster lookup

    # Debugging: Print number of unique headers read
    print(f"Number of unique headers read: {len(headers_to_keep)}")

    # Filter and write sequences
    with open(output_file, 'w') as filtered:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in headers_to_keep:
                filtered.write(f">{record.id}\n{record.seq}\n")
            else:
                # Debugging: Print IDs that are not found
                print(f"ID removed: {record.id}")

    # Check if all headers were used
    used_headers = set()
    with open(output_file, 'r') as check:
        for line in check:
            if line.startswith('>'):
                used_headers.add(line[1:].strip())

    # Debugging: Print headers not used
    unused_headers = headers_to_keep - used_headers
    print(f"Headers not used: {unused_headers}")

def main():
    parser = argparse.ArgumentParser(description="Extract contigs from a FASTA file based on headers from an annotation file.")
    parser.add_argument("--fasta_file", help="Input FASTA file")
    parser.add_argument("--annotation_file", help="Input annotation file")
    parser.add_argument("--output_file", help="Output FASTA file with filtered contigs")
    
    args = parser.parse_args()

    extractContigs(args.fasta_file, args.annotation_file, args.output_file)

if __name__ == "__main__":
    main()