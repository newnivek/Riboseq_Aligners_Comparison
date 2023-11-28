import sys

sam_file = sys.argv[1]
bed_file = sys.argv[2]

with open(sam_file, 'r') as sam, open(bed_file, 'w') as bed:
    for line in sam:
        if line.startswith('@'):
            continue  # Skip header lines

        fields = line.strip().split('\t')
        qname = fields[0]  # Query template NAME
        flag = int(fields[1])  # bitwise FLAG
        rname = fields[2]  # Reference sequence NAME
        pos = int(fields[3])  # 1-based leftmost mapping POSition
        mapq = fields[4]  # MAPping Quality

        # Skip unmapped reads
        if flag & 0x4:
            continue

        # Convert start position from 1-based to 0-based for BED format
        start = pos - 1
        # Estimate the end position (you might need to adjust this based on your data)
        end = start + len(fields[9])  # Using the length of the sequence

        bed.write(f"{rname}\t{start}\t{end}\t{qname}\t{mapq}\n")
