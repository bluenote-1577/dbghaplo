import pysam
import sys
import subprocess

# Define the input BAM file, output BAM file, and the IDs file
input_bam = sys.argv[1]
output_bam = input_bam + ".tagged.bam"
if len(sys.argv) > 2:
    ids_file = sys.argv[2]
else:
    ids_file = "dbghap_output/ids.txt"

print(ids_file)

# Read the ids.txt file and create a dictionary to map accessions to haplotypes
accession_to_hp = {}
with open(ids_file, 'r') as f:
    for line_number, line in enumerate(f):
        # Split the line into accessions and remove any whitespace
        accessions = line.strip().split()
        # Assign each accession an HP tag corresponding to the line number (1-indexed)
        contig = accessions[0]
        print(contig)
        hap = accessions[2].split(':')[1]
        ra = accessions[1].split(':')[1]
        tag = hap + '_' + ra
        for accession in accessions[3:]:
            accession_to_hp[accession] = tag

# Open the input BAM file for reading and the output BAM file for writing
with pysam.AlignmentFile(input_bam, "rb") as in_bam, pysam.AlignmentFile(output_bam, "wb", header=in_bam.header) as out_bam:
    for read in in_bam:
        # Extract the record accession from the read name
        accession = read.query_name
        # Check if the accession exists in the mapping
        if accession in accession_to_hp:
            if accession_to_hp[accession] != 'unassigned':
                # Add the HP tag to the read
                read.set_tag('HP', accession_to_hp[accession])
            # Write the read to the output BAM file
            out_bam.write(read)

# Index the output BAM file
subprocess.run(["samtools", "index", output_bam])
