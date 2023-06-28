aln_file = sys.argv[1]
bed_file = sys.argv[2]

with open(aln_file, 'r') as aln, open(bed_file, 'w') as bed:
    for line in aln:
        if line.startswith('>'):
            fields = line.strip().split('\t')
            identifier = fields[0][1:]  # remove the ">" at the start
            chromosome = identifier.split('-')[0]
            start = int(fields[2])   # convert to 0-based coordinate
            end = start + 28  # estimate end position based on read length
            bed.write(f"{chromosome}\t{start}\t{end}\n")
