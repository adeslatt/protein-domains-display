"""
####################################################################################################
# Script Name: generate_bed.py
# Description: Converts mapped protein domain coordinates into BED format for UCSC Genome Browser.
# Author: Anne Deslattes Mays, PhD
# Company: Science and Technology Consulting, LLC
# Copyright (C) 2025 Anne Deslattes Mays. All Rights Reserved.
# Date: 2025-02-16
####################################################################################################
"""

import pandas as pd
import os
import glob

def generate_bed(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    mapped_files = glob.glob(f"{input_dir}/*.tsv")

    for mapped_file in mapped_files:
        filename = os.path.basename(mapped_file).replace("_mapped_coordinates.tsv", "")
        output_file = f"{output_dir}/{filename}_domains.bed"

        df = pd.read_csv(mapped_file, sep="\t")
        df_bed = df[['chr', 'genomic_start', 'genomic_end', 'domain']]
        df_bed.to_csv(output_file, sep="\t", index=False, header=False)

        print(f"Saved BED: {output_file}")

if __name__ == "__main__":
    generate_bed("output/mapped_coordinates", "output")

