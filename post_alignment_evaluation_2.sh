#!/bin/bash

# Define variables for file paths
SAM_FILE="cds_aligned_hisat2_errFree.sam"  # Replace with the actual SAM file name
BAM_FILE="aligned.bam"
SORTED_BAM_FILE="aligned_sorted.bam"
BED_FILE="aligned.bed"
ORIGINAL_BED="original.bed"  # Replace with the actual BED file name
INTERSECT_BED="intersect.bed"

python3 sam2bed.py $SAM_FILE $BED_FILE

# Convert SAM to BAM
# samtools view -S -b $SAM_FILE > $BAM_FILE

# Sort the BAM file
# samtools sort $BAM_FILE -o $SORTED_BAM_FILE

# Index the sorted BAM file
# samtools index $SORTED_BAM_FILE

# Convert sorted BAM to BED
# bedtools bamtobed -i $SORTED_BAM_FILE > $BED_FILE

# Intersect the BED file with the original positions
bedtools intersect -F 1 -f 1 -a $BED_FILE -b $ORIGINAL_BED > $INTERSECT_BED


# Paths to the BED files for metrics calculation
aligned_bed=$BED_FILE
intersected_bed=$INTERSECT_BED

# Count the number of reads in each file
total_reads=$(wc -l < $ORIGINAL_BED)
aligned_reads=$(wc -l < $aligned_bed)
correctly_mapped=$(wc -l < $intersected_bed)

# Calculate TP, FP, FN
TP=$correctly_mapped
FP=$(expr $aligned_reads - $correctly_mapped)
FN=$(expr $total_reads - $correctly_mapped - $FP)

# Calculate Sensitivity (Recall), Precision, and F1 Score
Sensitivity=$(echo "$TP $FN" | awk '{print $1 / ($1 + $2)}')
Precision=$(echo "$TP $FP" | awk '{print $1 / ($1 + $2)}')
F1=$(echo "$Precision $Sensitivity" | awk '{print 2 * ($1 * $2) / ($1 + $2)}')

# Calculate Mapping Rate
MappingRate=$(echo "$aligned_reads $total_reads" | awk '{print $1 / $2}')

# Print the metrics
echo "Total reads: $total_reads"
echo "Aligned reads: $aligned_reads"
echo "Correctly Mapped (TP): $TP"
echo "Incorrectly Mapped (FP): $FP"
echo "Unmapped (FN): $FN"

echo "Sensitivity (Recall): $Sensitivity"
echo "Precision: $Precision"
echo "F1 Score: $F1"
echo "Mapping Rate: $MappingRate"
