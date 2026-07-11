from Bio import SeqIO
import pandas as pd
import sys
from joblib import Parallel, delayed
import os
import tempfile

# Parse input arguments
fasta_file = sys.argv[1]
dvf_file = sys.argv[2]
output_file = sys.argv[3]
temp_dirpath = sys.argv[4]

# Get DVF headers
df = pd.read_csv(dvf_file, sep='\t')
df = df[df['score'] > 0.95].dropna()
df['clean_name'] = df['name'].str.split().str[0]
headers_to_keep = set(df['clean_name'])

print(f"Number of unique headers to keep: {len(headers_to_keep)}")

# Function to filter FASTA records in a chunk
def process_chunk(chunk, headers_to_keep, output_chunk_file):
    with open(output_chunk_file, 'w') as filtered:
        for record in chunk:
            if record.id in headers_to_keep:
                filtered.write(f">{record.id}\n{record.seq}\n")

# Split FASTA into chunks
def split_fasta(fasta_file, num_chunks):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    chunk_size = max(1, len(records) // num_chunks)
    chunks = [records[i:i + chunk_size] for i in range(0, len(records), chunk_size)]
    return chunks

# Use custom temp directory
with tempfile.TemporaryDirectory(dir=temp_dirpath) as temp_dir:
    num_chunks = 20
    chunks = split_fasta(fasta_file, num_chunks)
    output_chunk_files = [os.path.join(temp_dir, f"output_chunk_{i}.fasta") for i in range(len(chunks))]

    # Process in parallel
    Parallel(n_jobs=num_chunks)(
        delayed(process_chunk)(chunk, headers_to_keep, out_file)
        for chunk, out_file in zip(chunks, output_chunk_files)
    )

    # Combine results
    with open(output_file, 'w') as outfile:
        for fname in output_chunk_files:
            with open(fname, 'r') as infile:
                outfile.write(infile.read())

# Check unused headers
used_headers = set()
with open(output_file, 'r') as check:
    for line in check:
        if line.startswith('>'):
            used_headers.add(line[1:].strip())

unused_headers = headers_to_keep - used_headers
print(f"Headers not used: {unused_headers}")
