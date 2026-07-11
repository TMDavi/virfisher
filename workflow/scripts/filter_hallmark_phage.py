import pandas as pd
from Bio import SeqIO
import os
import argparse

hallmark_genes = ["portal", "terminase", "major capsid protein"]

def findHallmark(annotation_description):
    """
    Check if the annotation contains any hallmark phage gene names.
    """
    if isinstance(annotation_description, str):
        lower_annotation = annotation_description.lower()
        for substring in hallmark_genes:
            if substring in lower_annotation:
                return True
    return False 

def countHallmarkGenes(annot_file):
    df = pd.read_csv(annot_file, sep='\t', encoding='utf-8-sig')
    df.columns = df.columns.str.strip()
    df['is_hallmark'] = df['annotation_description'].apply(findHallmark)
    df["contig_id"] = df['gene'].str.rsplit('_', n=1).str[0]

    df = df[['contig_id','is_hallmark']]
    df = df[df['is_hallmark'] == True].dropna()
    df = df.groupby('contig_id').size().reset_index(name='hallmark_count')
    df = df[df['hallmark_count'] > 0]
    
    return df

def filterHallmarkPhage(fasta_file, annot_file, out_file):
    hallmark_df = countHallmarkGenes(annot_file)
    hallmark_contigs = set(hallmark_df['contig_id'])

    records_to_keep = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in hallmark_contigs:
            records_to_keep.append(record)

    SeqIO.write(records_to_keep, out_file, "fasta")
    return hallmark_contigs  # Return list of contigs to use in other filters

def filterHallmarkPhageAnnotations(annot_file, contigs_to_keep, out_file):
    df = pd.read_csv(annot_file, sep='\t')
    df['contig_id'] = df['gene'].str.rsplit('_', n=1).str[0]
    df_filtered = df[df['contig_id'].isin(contigs_to_keep)]
    df_filtered.to_csv(out_file, sep='\t', index=False)
    return out_file

def filterHallmarkPhageQuality(qual_file, contigs_to_keep, out_file):
    quality_df = pd.read_csv(qual_file, sep="\t")
    quality_df['contig_id'] = quality_df['contig_id'].str.split(r"\|\|").str[0]
    filt_quality_df = quality_df[quality_df["contig_id"].isin(contigs_to_keep)]
    filt_quality_df.to_csv(out_file, sep='\t', index=False)
    return out_file

def main():
    parser = argparse.ArgumentParser(description="Filter FASTA, annotation, and quality files to keep only contigs with hallmark phage genes.")
    parser.add_argument("--fasta_file", required=True, help="Input FASTA file containing contigs.")
    parser.add_argument("--annot_file", required=True, help="Input annotation file in TSV format.")
    parser.add_argument("--qual_file", required=True, help="Input quality file in TSV format.")
    parser.add_argument("--out_dir", required=True, help="Output directory to save filtered files.")
    args = parser.parse_args()

    if not os.path.exists(args.fasta_file):
        raise FileNotFoundError(f"Input FASTA file '{args.fasta_file}' does not exist.")
    if not os.path.exists(args.annot_file):
        raise FileNotFoundError(f"Input annotation file '{args.annot_file}' does not exist.")
    if not os.path.exists(args.qual_file):
        raise FileNotFoundError(f"Input quality file '{args.qual_file}' does not exist.")
    
    os.makedirs(args.out_dir, exist_ok=True)

    # Filter and get contigs
    filtered_contigs = filterHallmarkPhage(
        args.fasta_file,
        args.annot_file,
        os.path.join(args.out_dir, "filtered_with_hallmark.fasta")
    )

    # Filter annotation
    filterHallmarkPhageAnnotations(
        args.annot_file,
        filtered_contigs,
        os.path.join(args.out_dir, "filtered_with_hallmark_annotations.tsv")
    )

    # Filter quality
    filterHallmarkPhageQuality(
        args.qual_file,
        filtered_contigs,
        os.path.join(args.out_dir, "filtered_with_hallmark_quality.tsv")
    )

if __name__ == "__main__":
    main()
