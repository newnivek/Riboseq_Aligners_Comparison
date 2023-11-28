import sys

def sam_to_fastq(sam_file, fastq_file):
    with open(sam_file, 'r') as sam, open(fastq_file, 'w') as fq:
        for line in sam:
            if line.startswith('@'):
                continue  # Skip header lines
            
            fields = line.strip().split('\t')
            qname = fields[0]  # Query template NAME
            seq = fields[9]    # segment SEQuence
            qual = fields[10]  # QUALity scores

            # Write to FastQ format
            fq.write(f'@{qname}\n{seq}\n+\n{qual}\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python sam_to_fastq.py <input.sam> <output.fq>")
        sys.exit(1)

    sam_file = sys.argv[1]
    fastq_file = sys.argv[2]
    sam_to_fastq(sam_file, fastq_file)
