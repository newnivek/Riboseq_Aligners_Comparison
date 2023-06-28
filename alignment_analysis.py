import pysam
from collections import defaultdict

def calculate_metrics(aligned_bam_path, intersected_bed_path):
    total_reads = 0
    correctly_mapped = 0
    incorrectly_other = 0
    unmapped = 0

    intersected_reads = set()
    with open(intersected_bed_path, 'r') as bed_file:
        for line in bed_file:
            fields = line.strip().split('\t')
            intersected_reads.add(fields[3])  # assuming read name is in 4th column

    with pysam.AlignmentFile(aligned_bam_path, 'rb') as aligned_file:
        for read in aligned_file:
            total_reads += 1
            if read.is_unmapped:
                unmapped += 1
            elif read.query_name in intersected_reads:
                correctly_mapped += 1
            else:
                incorrectly_other += 1

    precision = correctly_mapped / (correctly_mapped + incorrectly_other)
    recall = correctly_mapped / (correctly_mapped + incorrectly_other + unmapped)

    return {
        'Total reads': total_reads,
        'Correctly mapped': correctly_mapped,
        'Incorrectly mapped (other)': incorrectly_other,
        'Unmapped': unmapped,
        'Precision': precision,
        'Recall': recall,
    }

if __name__ == "__main__":
    aligned_bam_path = 'aligned_reads.bam'
    intersected_bed_path = 'intersected.bed'
    metrics = calculate_metrics(aligned_bam_path, intersected_bed_path)
    for metric, value in metrics.items():
        print(f'{metric}: {value}')
