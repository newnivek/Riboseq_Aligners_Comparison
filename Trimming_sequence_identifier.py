from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Input file
input_file = "GCF_003668045.3_CriGri-PICRH-1.0_cds_from_genomic.fna"

# Output file
output_file = "trimmed_identifiers.fna"

# This list will store your new records
new_records = []

# Open your fasta file
for record in SeqIO.parse(input_file, "fasta"):
    # Split the identifier at the "|" character and keep the second part (index 1)
    # Then split this part again at the "_cds" and keep the first part (index 0)
    new_id = record.id.split("|")[1].split("_cds")[0]
    
    # Create a new description removing the old identifier
    new_description = record.description.replace(record.id, "")
    
    # Create a new record with the new identifier, new description, and the old sequence
    new_record = SeqRecord(record.seq, id=new_id, description=new_description)
    
    # Add the new record to the list
    new_records.append(new_record)

# Write the new records to the output file
SeqIO.write(new_records, output_file, "fasta")
