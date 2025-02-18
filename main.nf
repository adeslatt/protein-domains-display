/*
####################################################################################################
# Pipeline Name: Protein Domain UCSC Visualization
# Description: Processes protein domain coordinates, maps to genome, and generates UCSC-compatible
#             bigBed and bigWig files for visualization.
# Author: Anne Deslattes Mays, PhD
# Company: Science and Technology Consulting, LLC
# Copyright (C) 2025 Anne Deslattes Mays. All Rights Reserved.
# Date: 2025-02-16
####################################################################################################
*/

nextflow.enable.dsl=2

// === Input Channels ===
gff_files = Channel.fromPath("${params.gff_dir}/*.gff").map { file(it) }
count_files = Channel.fromPath("${params.counts_dir}/*.csv").map { file(it) }

// === Process 1: Parse GFF Files ===
process parseGFF {
    input:
    path gff_files

    output:
    path("${params.output_dir}/parsed_gff/*.tsv")

    script:
    """
    python scripts/parse_gff.py ${params.gff_dir} ${params.output_dir}/parsed_gff
    """
}

// === Process 2: Process Protein Count Matrices ===
process processCounts {
    input:
    path count_files

    output:
    path("${params.output_dir}/processed_counts/*.tsv")

    script:
    """
    python scripts/process_counts.py ${params.counts_dir} ${params.output_dir}/processed_counts
    """
}

// === Process 3: Map AA to Genome ===
process mapAAToGenome {
    input:
    path gff_files
    path count_files

    output:
    path("${params.output_dir}/mapped_coordinates/*.tsv")

    script:
    """
    python scripts/map_aa_to_genome.py ${params.output_dir}/parsed_gff ${params.output_dir}/processed_counts ${params.output_dir}/mapped_coordinates
    """
}

// === Process 4: Generate BED Files ===
process generateBED {
    input:
    path("${params.output_dir}/mapped_coordinates/*.tsv")

    output:
    path("${params.output_dir}/*.bed")

    script:
    """
    python scripts/generate_bed.py ${params.output_dir}/mapped_coordinates ${params.output_dir}
    """
}

// === Process 5: Convert BED to bigBed ===
process convertToBigBed {
    input:
    path("${params.output_dir}/*.bed")

    output:
    path("${params.output_dir}/*.bb")

    script:
    """
    bedToBigBed ${params.output_dir}/*.bed ${params.chrom_sizes} ${params.output_dir}/*.bb
    """
}

// === Process 6: Generate bigWig ===
process generateBigWig {
    input:
    path("${params.output_dir}/mapped_coordinates/*.tsv")

    output:
    path("${params.output_dir}/*.bw")

    script:
    """
    bedGraphToBigWig ${params.output_dir}/mapped_coordinates/*.tsv ${params.chrom_sizes} ${params.output_dir}/*.bw
    """
}

// === Workflow Execution ===
workflow {
    gff_output = parseGFF(gff_files)
    counts_output = processCounts(count_files)
    mapped_coords = mapAAToGenome(gff_output, counts_output)
    bed_files = generateBED(mapped_coords)
    bigbed_files = convertToBigBed(bed_files)
    bigwig_files = generateBigWig(mapped_coords)
}

