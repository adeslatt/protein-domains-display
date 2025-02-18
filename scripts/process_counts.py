"""
####################################################################################################
# Script Name: process_counts.py
# Description: Reads and processes protein domain count matrices, renaming columns for consistency.
# Author: Anne Deslattes Mays, PhD
# Company: Science and Technology Consulting, LLC
# Copyright (C) 2025 Anne Deslattes Mays. All Rights Reserved.
# Date: 2025-02-16
####################################################################################################
"""

import pandas as pd
import glob
import os

def process_counts(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    csv_files = glob.glob(f"{input_dir}/*.csv")

    for csv_file in csv_files:
        protein_name = os.path.basename(csv_file).split("_")[0]
        output_file = f"{output_dir}/{protein_name}_processed_counts.tsv"

        df = pd.read_csv(csv_file)
        df.rename(columns={"Protein": "protein_id", "Domain": "domain", "Protein_AA_Coords": "aa_coords"}, inplace=True)
        df.to_csv(output_file, sep="\t", index=False)

        print(f"Processed {csv_file} -> {output_file}")

if __name__ == "__main__":
    process_counts("data/protein_matrices", "output/processed_counts")

