#!/usr/bin/env python3

from pathlib import Path
import argparse
import yaml

parser = argparse.ArgumentParser(
    description="Generate a configuration file for Virfisher."
)

parser.add_argument(
    "--input_folder",
    help="Directory containing paired FASTQ files."
)

parser.add_argument(
    "--output_folder",
    help="Directory where all results will be stored"
)

parser.add_argument(
    "--config",
    default="config.yaml",
    help="Output configuration file."
)

parser.add_argument(
    "--threads",
    type=int,
    default=4
)

parser.add_argument(
    "--mem_mb",
    type=int,
    default=256000
)

args = parser.parse_args()

input_dir = Path(args.input_folder)

if not input_dir.is_dir():
    raise SystemExit(f"{input_dir} is not a directory.")

config = {
    "WORKDIR":"workflow",
    "OUTDIR": args.output_folder,
    "resources": {
        "threads": args.threads,
        "mem_mb": args.mem_mb
    },
    "samples": {},
    "fastp": {
        "quality": 20
    },
    "metaspades": {
        "klist": "27,37,47,57,67,77,87,97,107,117,127",
        "mem_gb": 500
    }
}

for r1 in sorted(input_dir.rglob("*_R1.*")):

    sample = r1.name.split("_R1")[0]
    r2 = r1.with_name(r1.name.replace("_R1", "_R2"))

    if not r2.exists():
        print(f"Warning: missing pair for {r1.name}")
        continue

    config["samples"][sample] = {
        "forward": str(r1.resolve()),
        "reverseR": str(r2.resolve())
    }

with open(args.config, "w") as f:
    yaml.safe_dump(config, f, sort_keys=False)

print(f"Wrote {args.config}")