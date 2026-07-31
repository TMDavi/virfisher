import pandas as pd
import argparse 

def mergeDataFrame(filelist):
    df = pd.DataFrame()
    for file in filelist:
        sampleID = file.split('/')[-2]
        temp_df = pd.read_csv(file, sep='\t', names=['contig_id', 'length', 'Mapped', 'Unmapped'])
        temp_df = pd.DataFrame({
            'contig_id': temp_df['contig_id'],
            sampleID: temp_df['Mapped']
        })
        if df.empty:
            df = temp_df 
        else:
            df = pd.merge(df, temp_df, on='contig_id', how='outer')
    return df

def main():
    parser = argparse.ArgumentParser(description='Process mapping data files.')
    parser.add_argument('-i', '--input', nargs='+', help='List of sample files to process', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file name', default='idxstats_merged.tsv')
    args = parser.parse_args()

    df = mergeDataFrame(args.input)
    df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
