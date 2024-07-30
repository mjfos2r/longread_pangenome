#!/usr/local/bin/python3.11
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

def translate_nucleotide_to_protein(input_fasta, output_fasta):
    # Read the nucleotide sequences from the input FASTA file
    nucleotide_sequences = SeqIO.parse(input_fasta, "fasta")

    # Translate each nucleotide sequence to protein
    protein_sequences = []
    for nucleotide_sequence in tqdm(nucleotide_sequences, desc="sequences to translate"):
        protein_sequence = nucleotide_sequence.seq.translate(to_stop=True, table=11)
        protein_record = SeqRecord(protein_sequence, id=nucleotide_sequence.id, description="translated sequence")
        protein_sequences.append(protein_record)

    # Write the protein sequences to the output FASTA file
    SeqIO.write(protein_sequences, output_fasta, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python translate_fasta.py <input_fasta> <output_fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]

    translate_nucleotide_to_protein(input_fasta, output_fasta)
    print(f"Translation complete. Output written to {output_fasta}")
