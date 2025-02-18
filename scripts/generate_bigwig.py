"""
####################################################################################################
# Script Name: generate_bigwig.py
# Description: Converts protein domain count data to bigWig format for UCSC visualization.
#              Uses bedGraphToBigWig for the conversion.
# Author: Anne Deslattes Mays, PhD
# Company: Science and Technology Consulting, LLC
# Copyright (C) 2025 Anne Deslattes Mays. All Rights Reserved.
# Date: 2025-02-16
####################################################################################################
"""

import os
import glob
import pandas as pd
import subprocess

CHROM_SIZES = "data/hg38.chrom.sizes"

def convert_counts_to_bigwig(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    mapped_files = glob.glob(f"{input_dir}/*_mapped_coordinates.tsv")

    for mapped_file in mapped_files:
        filename = os.path.basename(mapped_file).replace("_mapped_coordinates.tsv", "")

        df = pd.read_csv(mapped_file, sep="\t")

        sample_columns = df.columns[5:]  # Assuming samples start at column index 5

        for sample in sample_columns:
            bedgraph_file = f"{output_dir}/{filename}_{sample}.bedGraph"
            bigwig_file = f"{output_dir}/{filename}_{sample}.bw"

            df_sample = df[['chr', 'genomic_start', 'genomic_end', sample]]
            df_sample.to_csv(bedgraph_file, sep="\t", index=False, header=False)

            cmd = f"bedGraphToBigWig {bedgraph_file} {CHROM_SIZES} {bigwig_file}"
            subprocess.run(cmd, shell=True, check=True)

            print(f"Converted {bedgraph_file} -> {bigwig_file}")

if __name__ == "__main__":
    convert_counts_to_bigwig("output/mapped_coordinates", "output")

