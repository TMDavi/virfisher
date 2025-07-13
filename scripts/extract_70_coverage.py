import pandas as pd

files = [
    'MP2_A_coverage.txt',
    'MP2_B_coverage.txt',
    'MP2_C_coverage.txt',
    'MP3_A_coverage.txt',
    'MP3_B_coverage.txt',
    'MP3_C_coverage.txt',
    'MP4_A_coverage.txt',
    'MP4_B_coverage.txt',
    'MP4_C_coverage.txt',
    'MP5_A_coverage.txt',
    'MP5_B_coverage.txt',
    'MP5_C_coverage.txt',
]

def get_70_coverage(file_path):
    """
    Parses the coverage file and returns a DataFrame with the relevant data.
    """
    df = pd.read_csv(file_path, sep='\t',header=None, skiprows=1, names=['contig','coverage'])
    df = df[df['coverage'] > 0.7]  # Filter for coverage greater than 70
    return df

def parse_coverage(file_list):
    df = pd.DataFrame()
    temp_df= pd.DataFrame()
    for file in file_list:
        temp_df = get_70_coverage(file)
        temp_df = temp_df.rename(columns={'contig': 'contig_id', 'coverage': f"{file.split('_coverage')[0]}"})
        if df.empty:
            df = temp_df 
        else:
            df = pd.merge(df, temp_df, on='contig_id', how='outer')
    
    df = df.fillna(0)  # Fill NaN values with 0
    return df

final= parse_coverage(files)
final.to_csv('contigs_coverage_morethan70.txt', sep='\t', index=False)