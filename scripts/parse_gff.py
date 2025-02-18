"""
####################################################################################################
# Script Name: parse_gff.py
# Description: Parses GFF files to extract protein-to-genome mapping and exon coordinates.
# Author: Anne Deslattes Mays, PhD
# Company: Science and Technology Consulting, LLC
# Copyright (C) 2025 Anne Deslattes Mays. All Rights Reserved.
# Date: 2025-02-16
####################################################################################################
"""

import pandas as pd
import glob
import os

def parse_gff(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    gff_files = glob.glob(f"{input_dir}/*.gff")

    for gff_file in gff_files:
        filename = os.path.basename(gff_file).replace(".gff", "")
        output_file = f"{output_dir}/{filename}_parsed_gff.tsv"

        gff_data = []
        with open(gff_file, "r") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if fields[2] == "CDS":
                    attr = {k: v for k, v in (item.split("=") for item in fields[8].split(";") if "=" in item)}
                    gff_data.append({
                        "chr": fields[0],
                        "start": int(fields[3]),
                        "end": int(fields[4]),
                        "strand": fields[6],
                        "protein_id": attr.get("protein_id", None),
                        "transcript_id": attr.get("Parent", None)
                    })

        df = pd.DataFrame(gff_data)
        df.to_csv(output_file, sep="\t", index=False)
        print(f"Parsed {gff_file} -> {output_file}")

if __name__ == "__main__":
    parse_gff("data/protein_domain_coordinates", "output/parsed_gff")

