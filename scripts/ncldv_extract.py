import pandas as pd
import argparse

def split_if_semicolon(s):
    if ";" in s:
        return s.split(";")[0]
    else:
        return s

def filter_by_length(string, min_length=1000):
    length = string.split("_")[3]
    if length.isdigit() and int(length) >= min_length:
        return string
    
def main():
    parser = argparse.ArgumentParser(description='Extract potential NCLDV contigs from a full output file.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input file path')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file path')
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t')
    potential = df['protein_ids'].apply(lambda x: split_if_semicolon(x))
    potential = potential.apply(lambda x: filter_by_length(x) if isinstance(x, str) else x)
    potential = potential.to_frame()
    potential = potential.dropna()
    potential = potential.rename(columns={'protein_ids': 'contig_id'})
    potential.to_csv(args.output, sep='\t', index=False)
