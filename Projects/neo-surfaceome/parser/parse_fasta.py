#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_fasta.py - A script to parse a multi-line FASTA file containing coding sequences into a TSV format.
"""

import sys

def parse_fasta(fasta_file, type, output_tsv):
    """
    Parses a multi-line FASTA file and extracts transcript ID, sequence length, and sequence.
    """
    try:
        with open(fasta_file, 'r') as infile:
            lines = infile.readlines()
            if not lines or not lines[0].startswith('>'):
                print("Error: Invalid FASTA format.")
                sys.exit(1)
    except FileNotFoundError:
        print(f"Error: File '{fasta_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    with open(fasta_file, 'r') as infile, open(output_tsv, 'w') as outfile:
        if type == "protein":
            outfile.write("protein_ID\tlength\tsequence\n")
        elif type == "transcript":
            outfile.write("transcript_ID\tlength\tsequence\n")
        else:
            print("Error: Type must be either 'protein' or 'transcript'.")
            sys.exit(1)
        
        _id = None
        sequence = []
        
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if _id and sequence:
                    seq = ''.join(sequence)
                    outfile.write(f"{_id}\t{len(seq)}\t{seq}\n")
                
                parts = line[1:].split()
                _id = parts[0] if parts else None
                sequence = []
            else:
                sequence.append(line)
        
        if _id and sequence:  # Write last entry
            seq = ''.join(sequence)
            outfile.write(f"{_id}\t{len(seq)}\t{seq}\n")
    
    print(f"Parsing complete. Output saved to {output_tsv}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python parse_fasta.py <input.fasta> <type='protein'/'transcript'> <output.tsv>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    type = sys.argv[2]
    output_tsv = sys.argv[3]
    parse_fasta(input_fasta, type, output_tsv)
