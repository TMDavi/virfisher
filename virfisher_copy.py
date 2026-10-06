#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

#RUN MODE

parser = argparse.ArgumentParser(prog="virfisher")
subparsers = parser.add_subparsers(dest="command", required=True)

run_parser = subparsers.add_parser(
    "run",
    help="Run the viral metagenomics pipeline"
)

run_parser.add_argument(
    "--input_folder",
    help="Directory containing paired FASTQ files."
)

run_parser.add_argument(
    "--input_sheet",
    help="Sheet containg fastq file paths"
)

run_parser.add_argument(
    "--outdir",
    default="results",
    help="Output directory."
)

run_parser.add_argument(
    "--threads",
    type=int,
    default=4,
    help="Maximum number of threads."
)

run_parser.add_argument(
    "--mem_mb",
    type=int,
    default=256000,
    help="Maximum memory (MB) available per job."
)

run_parser.add_argument(
    "--config",
    default="config.yaml",
    help="Configuration file to generate."
)

run_parser.add_argument(
    "--dry_run",
    default=False,
    help="Check if the pipeline is generating all input files"
)


args_run = run_parser.parse_args()

# Step 1: Generate config.yaml

print("Generating configuration file...")

create_config = SCRIPT_DIR / "create-config.py"

if args_run.input_folder:

    create_config_cmd = [
        sys.executable,
        str(create_config),
        "--input_folder", args_run.input_folder,
        "--output_folder", args_run.outdir,
        "--config", args_run.config,
        "--threads", str(args_run.threads),
        "--mem_mb", str(args_run.mem_mb),
    ]

elif args_run.input_sheet:
    create_config_cmd = [
            sys.executable,
            str(create_config),
            "--input_sheet", args_run.input_sheet,
            "--output_folder", args_run.outdir,
            "--config", args_run.config,
            "--threads", str(args_run.threads),
            "--mem_mb", str(args_run.mem_mb),
        ]

subprocess.run(create_config_cmd, check=True)

# Step 2: Create output directory

Path(args_run.outdir).mkdir(parents=True, exist_ok=True)

# Step 3: Run Snakemake

print("Running Snakemake...")

snakemake_cmd = [
    "snakemake",
    "--snakefile", "workflow/main.smk",
    "--configfile", args_run.config,
    "--config", f"outdir={args_run.outdir}",
    "--cores", str(args_run.threads),
    "--use-conda"
]

snakemake_cmd_dryrun = [
    "snakemake",
    "--snakefile", "workflow/main.smk",
    "--configfile", args_run.config,
    "--config", f"outdir={args_run.outdir}",
    "--cores", str(args_run.threads),
    "--use-conda",
    "--dry-run"
]

if args_run.dry_run:
    subprocess.run(snakemake_cmd_dryrun, check=True)
else:
    subprocess.run(snakemake_cmd, check=True)

print("Pipeline completed successfully.")

#SETUP DATABASES

db_parser = subparsers.add_parser(
    "download_db",
    help="Run the viral metagenomics pipeline"
)
db_parser.add_argument(
    "--db_dir",
    help="Directory containing paired FASTQ files."
)

snakemake_cmd_dryrun = [
    "snakemake",
    "--snakefile", "workflow/rules/database_setup/setup_DBs.smk",
    "--cores", str(args_run.threads),
    "--use-conda"
]
