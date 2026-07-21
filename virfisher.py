#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

parser = argparse.ArgumentParser(
    description="Run the Virfisher workflow."
)

parser.add_argument(
    "--input",
    required=True,
    help="Directory containing paired FASTQ files."
)

parser.add_argument(
    "--outdir",
    default="results",
    help="Output directory."
)

parser.add_argument(
    "--threads",
    type=int,
    default=4,
    help="Maximum number of threads."
)

parser.add_argument(
    "--mem_mb",
    type=int,
    default=256000,
    help="Maximum memory (MB) available per job."
)

parser.add_argument(
    "--config",
    default="config.yaml",
    help="Configuration file to generate."
)

parser.add_argument(
    "--dry_run",
    default=False,
    help="Check if the pipeline is generating all input files"
)


args = parser.parse_args()

# Step 1: Generate config.yaml

print("Generating configuration file...")

create_config = SCRIPT_DIR / "create-config.py"

create_config_cmd = [
    sys.executable,
    str(create_config),
    "--input_folder", args.input,
    "--output_folder", args.outdir,
    "--config", args.config,
    "--threads", str(args.threads),
    "--mem_mb", str(args.mem_mb),
]

subprocess.run(create_config_cmd, check=True)

# Step 2: Create output directory

Path(args.outdir).mkdir(parents=True, exist_ok=True)

# Step 3: Run Snakemake

print("Running Snakemake...")

snakemake_cmd = [
    "snakemake",
    "--snakefile", "workflow/main.smk",
    "--configfile", args.config,
    "--config", f"outdir={args.outdir}",
    "--cores", str(args.threads),
    "--use-conda"
]

snakemake_cmd_dryrun = [
    "snakemake",
    "--snakefile", "workflow/main.smk",
    "--configfile", args.config,
    "--config", f"outdir={args.outdir}",
    "--cores", str(args.threads),
    "--use-conda",
    "--dry-run"
]

if args.dry_run:
    subprocess.run(snakemake_cmd_dryrun, check=True)
else:
    subprocess.run(snakemake_cmd, check=True)

print("Pipeline completed successfully.")