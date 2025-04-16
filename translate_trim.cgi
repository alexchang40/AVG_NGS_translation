#!/usr/bin/python3

import io
import cgi
import cgitb
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

# Enable CGI error reporting
cgitb.enable()

# Function to trim DNA to a multiple of 3
def trim_to_multiple_of_three(dna_seq):
    remainder = len(dna_seq) % 3
    if remainder != 0:
        return dna_seq[:-remainder]
    return dna_seq

# Function to translate DNA to protein in all six reading frames
def translate_dna(dna_seq):
    translations = []
    for strand, nuc in [(1, dna_seq), (-1, dna_seq.reverse_complement())]:
        nuc = trim_to_multiple_of_three(nuc)
        for frame in range(3):
            protein = nuc[frame:].translate(to_stop=False)
            translations.append((strand, frame, protein))
    return translations

# Function to trim the protein sequence
def trim_sequence(protein_seq, start_pattern, stop_pattern):
    start_match = start_pattern.search(protein_seq)
    stop_match = stop_pattern.search(protein_seq)

    if start_match and stop_match:
        start_index = start_match.start()
        stop_index = stop_match.end()
        if start_index < stop_index:
            return protein_seq[start_index:stop_index]
    return None

# Main CGI logic
print("Content-Type: text/html\n")  # HTTP header

# Parse form data
form = cgi.FieldStorage()

# Get uploaded file and other inputs
if "input_file" not in form or "output_file" not in form or "input_count" not in form or "start_pattern" not in form or "stop_pattern" not in form:
    print("<h1>Error: Missing form data</h1>")
    exit()

input_file = io.TextIOWrapper(form["input_file"].file, encoding = 'ascii')
output_file = form["output_file"].value
input_count = int(form["input_count"].value)

# Start and stop patterns
start_pattern = form["start_pattern"].value.upper()
stop_pattern = form["stop_pattern"].value.upper()

unique_proteins = {}
total_sequences = 0
total_matches = 0
unique_clones = 0

# Parse the FASTQ file
for record in SeqIO.parse(input_file, "fastq"):
    total_sequences += 1
    trimmed_dna = trim_to_multiple_of_three(record.seq)
    translations = translate_dna(trimmed_dna)

    for strand, frame, protein in translations:
        protein_seq = str(protein)
        trimmed_seq = trim_sequence(protein_seq, start_pattern, stop_pattern)

        if trimmed_seq:
            total_matches+=1    
            if trimmed_seq in unique_proteins:
                unique_proteins[trimmed_seq]["count"] += 1
            else:
                unique_proteins[trimmed_seq] = {
                    "count": 1,
                    "record": SeqRecord(
                        Seq(trimmed_seq),
                        id=record.id,
                        description=f"strand = {strand}, frame = {frame}"
                    )
                }
            break

unique_clones = len(unique_proteins)
clones_above_threshold = {
    seq: data for seq, data in unique_proteins.items()
    if data["count"] >= input_count
}

num_clones_above_thershold = len(clones_above_threshold)
total_counts_above_threshold = sum(data["count"] for data in clones_above_threshold.values())

# Write the output FASTA file
output_path = f"/var/www/html/{output_file}"  # Save output in web-accessible directory
with open(output_path, "w") as output_handle:
    for seq, data in unique_proteins.items():
        if data["count"] > (input_count - 1):
            data["record"].description += f", count = {data['count']}"
            SeqIO.write(data["record"], output_handle, "fasta-2line")

# Display results
print("<h1>Translation and Trimming Complete</h1>")
print(f"<p>Output file: <a href='/{output_file}'>{output_file}</a></p>")
print("<h2>Results</h2>")
if input_count == 1:
    print("<p>Threshold: all counts</p>")
else:
    print("<p>Threshold: {input_count}</p>")
print(f"""
<p>Total sequences processed: {total_sequences}</p>
<p>Total unique clones found: {unique_clones}</p>
<p>Unique clones above threshold: {num_clones_above_threshold}</p>
<p>Total count of threshold-passing clones: {total_counts_above_threshold}</p>
""")
