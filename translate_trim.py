#!/usr/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
import re

#start pattern is "MQLL", stop pattern is "VTVSS"
start_pattern = re.compile("MQLL")
stop_pattern = re.compile("VTVSS")

#Prompt the user for input and output file names
input_file = input("Enter the input file to be read: ")
output_file = input("Enter the name of the output fasta file: ")

#ensure that the DNA is a multiple of 3
def trim_to_multiple_of_three(dna_seq):
    remainder = len(dna_seq)%3
    if remainder != 0:
        return dna_seq[:-remainder]
    return dna_seq

#translate DNA to protein in all six reading frames
def translate_dna(dna_seq):
    #empty list of the translations
    translations = []
    for strand, nuc in [(1, dna_seq), (-1, dna_seq.reverse_complement())]:
        nuc = trim_to_multiple_of_three(nuc)
        for frame in range(3):
            protein = nuc[frame:].translate(to_stop=False)
            translations.append((strand, frame, protein))
    return translations

#trim the aa sequence 
def trim_sequence(protein_seq, start_pattern, stop_pattern):
    start_match = start_pattern.search(protein_seq)
    stop_match = stop_pattern.search(protein_seq)

    if start_match and stop_match:
        start_index = start_match.start()
        stop_index = stop_match.end()
        if start_index<stop_index:
            return protein_seq[start_index:stop_index]
    return None

unique_proteins = {}

#parse through the fastq
for record in SeqIO.parse(input_file, "fastq"):
    #debug
#    print(f"Original sequence: {record.seq}")
    
    trimmed_dna = trim_to_multiple_of_three(record.seq)
    #debug
#    print(f"Trimmed sequence: {trimmed_dna}, Length: {len(trimmed_dna)}")
    
    translations = translate_dna(trimmed_dna)
    for strand, frame, protein in translations:
        #debug
#        print(f"Translated protein (strand = {strand}, frame = {frame}): {protein}")

        protein_seq = str(protein)
        trimmed_seq = trim_sequence(protein_seq, start_pattern, stop_pattern)

        if trimmed_seq:
            if trimmed_seq in unique_proteins:
                unique_proteins[trimmed_seq]["count"] += 1
            else:
                unique_proteins[trimmed_seq] = {
                    "count": 1,
                    "record": SeqRecord(
                        Seq(trimmed_seq),
                        id = record.id,
                        description = f"strand = {strand}, frame = {frame}"
                    )
                }

totalcount = 0
#writes to the output file
with open(output_file, "w") as output_handle:
    for seq, data in unique_proteins.items():
        if data["count"]>0: #only write sequences with count >5    
            data["record"].description += f", count = {data['count']}"
            SeqIO.write(data["record"], output_handle, "fasta")
            totalcount += data['count']

print(f"Total count is {totalcount}")



