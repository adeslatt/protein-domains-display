"""
####################################################################################################
# Script Name: convert_to_bigbed.py
# Description: Converts BED files to bigBed format for UCSC visualization using bedToBigBed.
# Author: Anne Deslattes Mays, PhD
# Company: Science and Technology Consulting, LLC
# Copyright (C) 2025 Anne Deslattes Mays. All Rights Reserved.
# Date: 2025-02-16
####################################################################################################
"""

import os
import glob
import subprocess

CHROM_SIZES = "data/hg38.chrom.sizes"

def convert_bed_to_bigbed(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    bed_files = glob.glob(f"{input_dir}/*_domains.bed")

    for bed_file in bed_files:
        filename = os.path.basename(bed_file).replace("_domains.bed", "")
        output_file = f"{output_dir}/{filename}_domains.bb"

        cmd = f"bedToBigBed {bed_file} {CHROM_SIZES} {output_file}"
        subprocess.run(cmd, shell=True, check=True)
        print(f"Converted {bed_file} -> {output_file}")

if __name__ == "__main__":
    convert_bed_to_bigbed("output", "output")

