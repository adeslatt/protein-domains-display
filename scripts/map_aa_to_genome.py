"""
####################################################################################################
# Script Name: map_aa_to_genome.py
# Description: Maps protein amino acid coordinates to genomic coordinates using parsed GFF files.
# Author: Anne Deslattes Mays, PhD
# Company: Science and Technology Consulting, LLC
# Copyright (C) 2025 Anne Deslattes Mays. All Rights Reserved.
# Date: 2025-02-16
####################################################################################################
"""

import pandas as pd
import re
import glob
import os

def map_aa_to_genome(gff_dir, counts_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    gff_files = glob.glob(f"{gff_dir}/*.tsv")
    count_files = glob.glob(f"{counts_dir}/*.tsv")

    for gff_file in gff_files:
        filename = os.path.basename(gff_file).replace("_parsed_gff.tsv", "")
        count_file = f"{counts_dir}/{filename}_processed_counts.tsv"

        if not os.path.exists(count_file):
            print(f"Skipping {filename}: No matching counts file found.")
            continue

        gff_df = pd.read_csv(gff_file, sep="\t")
        df = pd.read_csv(count_file, sep="\t")

        def map_aa(protein_id, aa_coords):
            protein_exons = gff_df[gff_df["protein_id"] == protein_id].sort_values(["transcript_id", "start"])
            if protein_exons.empty:
                return None, None, None

            match = re.match(r"(\d+)-(\d+)", str(aa_coords))
            if not match:
                return None, None, None
            start_aa, end_aa = map(int, match.groups())

            best_match = None
            for transcript_id in protein_exons["transcript_id"].unique():
                transcript_exons = protein_exons[protein_exons["transcript_id"] == transcript_id]
                
                aa_pos = 1
                genomic_start, genomic_end = None, None

                for _, exon in transcript_exons.iterrows():
                    exon_length = (exon["end"] - exon["start"] + 1) // 3

                    if aa_pos <= start_aa < aa_pos + exon_length:
                        offset = (start_aa - aa_pos) * 3
                        genomic_start = exon["start"] + offset

                    if aa_pos <= end_aa < aa_pos + exon_length:
                        offset = (end_aa - aa_pos) * 3
                        genomic_end = exon["start"] + offset
                        break

                    aa_pos += exon_length
                
                if genomic_start and genomic_end:
                    best_match = (transcript_id, genomic_start, genomic_end)
                    break

            return best_match if best_match else (None, None, None)

        df[['transcript_id', 'genomic_start', 'genomic_end']] = df.apply(
            lambda row: map_aa(row['protein_id'], row['aa_coords']), axis=1, result_type="expand"
        )

        output_file = f"{output_dir}/{filename}_mapped_coordinates.tsv"
        df.to_csv(output_file, sep="\t", index=False)
        print(f"Mapped {filename} -> {output_file}")

if __name__ == "__main__":
    map_aa_to_genome("output/parsed_gff", "output/processed_counts", "output/mapped_coordinates")

