from Bio import SeqIO
import sys
from joblib import Parallel, delayed
import os
import tempfile
import pandas as pd
import argparse
import shutil 

def process_chunk(chunk, headers_to_keep, output_chunk_file):
    try:
        with open(output_chunk_file, 'w') as filtered:
            for record in chunk:
                if record.id in headers_to_keep:
                    
                    SeqIO.write(record, filtered, "fasta")
    except IOError as e:
        print(f"Error writing to chunk file {output_chunk_file}: {e}", file=sys.stderr)
        raise 

def split_fasta(fasta_file, num_chunks):
    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))
        if not records:
            print(f"Warning: No records found in {fasta_file}.", file=sys.stderr)
            return []
        chunk_size = max(1, len(records) // num_chunks)
        chunks = [records[i:i + chunk_size] for i in range(0, len(records), chunk_size)]
        return chunks
    except FileNotFoundError:
        print(f"Error: Input FASTA file not found at {fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing FASTA file {fasta_file}: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Filter FASTA records based on a list of headers using parallel processing.")
    parser.add_argument("-i","--fasta_file",help="Path to the input FASTA file.")
    parser.add_argument("-hd","--headers_file",help="Path to a file containing headers to keep. Only the first column will be used.")
    parser.add_argument("-o", "--output_file",help="Path to the output filtered FASTA file.")
    parser.add_argument("--temp_dir_parent",default=".",help="Parent directory for creating a unique temporary directory. ""A new, uniquely named temporary directory will be created inside this path. " "(default: current working directory)")
    parser.add_argument("-c","--num_chunks",type=int,default=4,help="Number of chunks to split the FASTA file into for parallel processing. ""Increasing this may increase memory usage if the FASTA is large. (default: 20)")
    parser.add_argument("-t","--n_jobs",type=int,default=-1,help="Number of parallel jobs to run. -1 means use all available CPU cores. (default: -1)")

    args = parser.parse_args()

    # --- Input Validation and Setup ---
    if not os.path.exists(args.fasta_file):
        print(f"Error: Input FASTA file not found at '{args.fasta_file}'", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(args.headers_file):
        print(f"Error: Headers file not found at '{args.headers_file}'", file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(args.temp_dir_parent):
        try:
            os.makedirs(args.temp_dir_parent)
            print(f"Created temporary directory parent: {args.temp_dir_parent}")
        except OSError as e:
            print(f"Error: Could not create temporary directory parent '{args.temp_dir_parent}': {e}", file=sys.stderr)
            sys.exit(1)

    # Read headers to filter
    headers = pd.read_csv(args.headers_file, sep='\t')
    headers = headers.iloc[:, 0].tolist() 
    headers_to_keep = set()
    try:
        for header in headers:
            if header: # Ensure not to add empty strings
                headers_to_keep.add(header)
    except IOError as e:
        print(f"Error reading headers file '{args.headers_file}': {e}", file=sys.stderr)
        sys.exit(1)

    if not headers_to_keep:
        print("Warning: No headers found in the headers file. Output FASTA will be empty.", file=sys.stderr)

    print(f"Number of unique headers to keep: {len(headers_to_keep)}")

    # Use a custom temporary directory (a unique one will be created inside temp_dir_parent)
    temp_dir_obj = None
    try:
        temp_dir_obj = tempfile.TemporaryDirectory(dir=args.temp_dir_parent)
        temp_dir = temp_dir_obj.name
        print(f"Using temporary directory: {temp_dir}")

        chunks = split_fasta(args.fasta_file, args.num_chunks)
        if not chunks:
            print("No FASTA records to process. Exiting.", file=sys.stderr)
            return

        # Generate unique temporary file names for each chunk's output
        output_chunk_files = [os.path.join(temp_dir, f"output_chunk_{i}.fasta") for i in range(len(chunks))]

        # Process in parallel
        print(f"Processing {len(chunks)} chunks using {args.n_jobs} parallel jobs...")
        Parallel(n_jobs=args.n_jobs)(
            delayed(process_chunk)(chunk, headers_to_keep, out_file)
            for chunk, out_file in zip(chunks, output_chunk_files)
        )
        print("Parallel processing complete. Combining results...")

        # Combine results
        with open(args.output_file, 'w') as outfile:
            for fname in output_chunk_files:
                if os.path.exists(fname): # Ensure chunk file exists before trying to read
                    with open(fname, 'r') as infile:
                        # Use shutil.copyfileobj for efficient file copying
                        shutil.copyfileobj(infile, outfile)
                else:
                    print(f"Warning: Temporary chunk file '{fname}' not found. It might have been empty or an error occurred.", file=sys.stderr)
        print(f"Filtered FASTA saved to: {args.output_file}")

        # Check unused headers (headers that were requested but not found/written)
        used_headers = set()
        if os.path.exists(args.output_file):
            try:
                with open(args.output_file, 'r') as check:
                    for record in SeqIO.parse(check, "fasta"): # Use SeqIO.parse for robustness
                        used_headers.add(record.id)
            except IOError as e:
                print(f"Error reading output file for unused headers check: {e}", file=sys.stderr)
        else:
            print(f"Warning: Output file '{args.output_file}' not found for checking used headers.", file=sys.stderr)

        unused_headers = headers_to_keep - used_headers
        if unused_headers:
            print(f"\n--- Headers from '{args.headers_file}' NOT found in the input FASTA or not written to output ({len(unused_headers)}): ---")
            for header in sorted(list(unused_headers)):
                print(f"- {header}")
            print("------------------------------------------------------------------------------------------------")
        else:
            print("\nAll specified headers were found and written to the output file.")

    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        # Clean up the temporary directory
        if temp_dir_obj:
            try:
                temp_dir_obj.cleanup()
                print(f"Cleaned up temporary directory: {temp_dir_obj.name}")
            except Exception as e:
                print(f"Error cleaning up temporary directory '{temp_dir_obj.name}': {e}", file=sys.stderr)

if __name__ == "__main__":
    main()