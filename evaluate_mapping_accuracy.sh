#!/bin/bash

# Input files
aligned_bam=$1
coding_regions_bed=$2
mapQ_threshold=$3

# Filter the aligned reads by the mapQ threshold
samtools view -b -q $mapQ_threshold $aligned_bam > output_sorted_filtered.bam

# Convert the filtered BAM file to BED format
bedtools bamtobed -i output_sorted_filtered.bam > aligned_reads_filtered.bed

# Find the intersecting reads between the filtered reads and the coding regions
bedtools intersect -a aligned_reads_filtered.bed -b $coding_regions_bed -wa -wb > intersected_filtered_output.bed

# Count the number of correctly mapped reads (intersecting), total reads, and total filtered reads
correctly_mapped=$(wc -l < intersected_filtered_output.bed)
total_reads=$(samtools view -c $aligned_bam)
total_filtered_reads=$(samtools view -c output_sorted_filtered.bam)

# Calculate incorrectly mapped and unmapped reads
incorrectly_mapped=$((total_filtered_reads - correctly_mapped))
unmapped=$((total_reads - total_filtered_reads))

# Calculate precision, recall, and F0.25-measure
precision=$(echo "scale=4; $correctly_mapped / ($correctly_mapped + $incorrectly_mapped)" | bc)
recall=$(echo "scale=4; $correctly_mapped / ($correctly_mapped + $unmapped)" | bc)
F_measure=$(echo "scale=4; (1.25 * $precision * $recall) / (0.25 * $precision + $recall)" | bc)

# Print results
echo "mapQ threshold: $mapQ_threshold"
echo "Precision: $precision"
echo "Recall: $recall"
echo "F0.25 measure: $F_measure"
echo
