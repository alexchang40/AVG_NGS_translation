#!/usr/bin/python3

import cgi
import json
import os
from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
import re



#Function to ensure that the DNA is a multiple of 3
def trim_to_multiple_of_three(dna_seq):
    remainder = len(dna_seq)%3
    if remainder != 0:
        return dna_seq[:-remainder]
    return dna_seq

#Function to translate DNA to protein in all six reading frames
def translate_dna(dna_seq):
    #empty list of the translations
    translations = []
    for strand, nuc in [(1, dna_seq), (-1, dna_seq.reverse_complement())]:
        nuc = trim_to_multiple_of_three(nuc)
        for frame in range(3):
            protein = nuc[frame:].translate(to_stop=False)
            translations.append((strand, frame, protein))
    return translations

#Function to trim the aa sequence based on start and stop patterns
def trim_sequence(protein_seq, start_pattern, stop_pattern):
    start_match = start_pattern.search(protein_seq)
    stop_match = stop_pattern.search(protein_seq)

    if start_match and stop_match:
        start_index = start_match.start()
        stop_index = stop_match.end()
        if start_index<stop_index:
            return protein_seq[start_index:stop_index]
    return None


#Main CGI logic
def main():
    print("Content-Type: application/json\n")

    form = cgi.FieldStorage()

    #check if the form was submitted
    if "file" in form and "start_seq" in form and "stop_seq" in form and "min_count" in form:
        #Get user inputs
        start_seq = form.getvalue("start_seq")
        stop_seq = form.getvalue("stop_seq")
        min_count = int(form.getvalue("min_count"))
        file_item = form["file"]

        #Save the uploaded file
        input_file = "uploaded_file.fastq"
        with open(input_file, "wb") as f:
            f.write(file_item.file.read())
        
        #compile regex patterns
        start_pattern = re.compile(start_seq)
        stop_pattern = re.compile(stop_seq)

        unique_proteins = {}
        total_sequences = 0   #counter for total sequences parsed
        unique_clones = 0     #counter for unique clones
 
        #process the FASTQ file
        for record in SeqIO.parse(input_file, "fastq"):
            total_sequences +=1 #increment total sequences counter
            translations = translate_dna(record.seq)
            for strand, frame, protein in translations:
                protein_seq = str(protein)
                trimmed_seq = trim_sequence(protein_seq, start_pattern, stop_pattern)
                
                if trimmed_seq:
                    if trimmed_seq in unique_proteins:
                        unique_proteins[trimmed_seq]["count"]+=1
                    else:
                        unique_clones += 1 #increment unique clones counter
                        unique_proteins[trimmed_seq] = {
                            "count": 1,
                            "record": SeqRecord(
                                Seq(trimmed_seq),
                                id = record.id,
                                description = f"strand = {strand}, frame = {frame}"
                            )
                        }
        #write output to FASTA file
        output_file = "output.fasta"
        total_count = 0
        with open(output_file, "w") as output_handle:
            for seq, data in unique_proteins.items():
                if data["count"] >= min_count:
                    data["record"].description += f", count = {data['count']}"
                    SeqIO.write(data["record"], output_handle, "fasta")
                    total_count += data['count']

            return render_template(
                "results.html",
                total_sequences = total_sequences,
                unique_clones = unique_clones,
                total_count = total_count,
                output_file = output_file
            )
        return render_template("index.html")

@app.route("/download<filename>")
def download(filename):
    return send_file(filename, as_attachment=True)

if __name__ == "__main__":
    main()
