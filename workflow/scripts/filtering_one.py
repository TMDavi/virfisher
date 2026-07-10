from Bio import SeqIO
import os
import pandas as pd
import argparse
import tempfile

def generateFilterlist(quality_file):
   
    #Step 1 - Apply filters to the quality summary file
    #get contig >=3000 and viral genes > 0
    df = pd.read_csv(quality_file,sep='\t')
    df1 = df[(df['contig_length'] >= 3000) & (df['viral_genes'] > 0)]
    print(len(df1))
    #get contig >=5000 viral genes ==0 and host genes = 0
    df = pd.read_csv(quality_file,sep='\t')
    df2 = df[(df['contig_length'] > 5000) & (df['viral_genes'] == 0) & (df['host_genes'] == 0)]
    print(len(df2))
    #get contigs >=10000 and viral genes == 0 and host genes = 1
    df = pd.read_csv(quality_file,sep='\t')
    df3 = df[(df['contig_length'] > 10000) & (df['viral_genes'] == 0) & (df['host_genes'] == 1)]
    print(len(df3))
    #Combine dataframes
    frames = [df1, df2, df3]
    df = pd.concat(frames)

    #Step 2 - Remove contigs with 3 times more host genes than viral genes and no proviral signals
    df_result = df.drop(df[(df['viral_genes'] * 3 < df['host_genes']) & (df['provirus'] == 'No')].index)
    
    df_result['contig_id'] = df_result['contig_id'].str.split(r"\|\|").str[0]  # Extract the part before '||'

    result = df_result['contig_id'].tolist()
    

    return result
    #filter.to_csv(f'MP{s}_filter1_list.txt',header=False, index=False)pi


def extract_contigs(fasta_file, headers_list, output_file):

    # Read headers to filter
    headers_to_keep = headers_list

    # Debugging: Print number of unique headers read
    print(f"Number of unique headers read: {len(headers_to_keep)}")

    # Filter and write sequences
    with open(output_file, 'w') as filtered:
        for record in SeqIO.parse(fasta_file, "fasta"):
            
            if record.id in headers_to_keep:
                filtered.write(f">{record.id}\n{record.seq}\n")
            else:
                # Debugging: Print IDs that are not found
                print(f"Did not meet quality parameters: {record.id}")


    # Check if all headers were used
    headers_to_keep = set(headers_to_keep)
    used_headers = set()
    with open(output_file, 'r') as check:
        for line in check:
            if line.startswith('>'):
                used_headers.add(line[1:].strip().split()[0])

    # Debugging: Print headers not used
    unused_headers = headers_to_keep - used_headers
    print(f"{len(unused_headers)} Headers not used: {unused_headers}")
    

def getLength(fasta):
    
    contig_length = pd.DataFrame(
        {'contig': [record.id for record in SeqIO.parse(fasta, 'fasta')],
         'length': [len(record.seq) for record in SeqIO.parse(fasta, 'fasta')]
        }
    )

    return contig_length

def main():
    parser = argparse.ArgumentParser(description="Extract contigs from a FASTA file based on headers.")
    parser.add_argument("--fasta_file", help="Input FASTA file")
    parser.add_argument("--quality_file", help="Input file containing checkv quality summary")
    parser.add_argument("--tempdir", help="Temporary directory")
    parser.add_argument("--output_file", help="Output filtered FASTA file")
    args = parser.parse_args()

    #Step1 = FGenerate the filter list based on checkV results
    headers_list = generateFilterlist(args.quality_file)
    print(headers_list)

    with tempfile.NamedTemporaryFile(mode='w+', dir=args.tempdir, delete=False, suffix=".fa") as temp_fasta:
        temp_fasta_name = temp_fasta.name

        # Step 3: Filter the contigs
        extract_contigs(args.fasta_file, headers_list, temp_fasta_name)
    
    
    #Step3 = CheckV displays the contig length before the trimming so it is expected to stil have contigs < 3000
    #We now check the contig length of the filtered contigs
    c_len = getLength(temp_fasta_name)

    #Step 4 = Finally remove contigs with less than 3000 bp with will not be annotated by geNomad
    more_than_3000 = c_len[c_len['length'] >= 3000]
    more_than_3000 = more_than_3000['contig'].tolist()

    extract_contigs(temp_fasta_name,more_than_3000, args.output_file)

    os.remove(temp_fasta_name)  

    
if __name__ == "__main__":
    main()