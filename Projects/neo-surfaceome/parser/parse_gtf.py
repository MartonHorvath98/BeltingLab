#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parse GTF file and extract relevant information for the database annotation of each transcript.
"""

import sys
import re


def parse_gtf(gtf_file, output_tsv):
    """
    Parses the Ensembl GTF file and extracts transcript ID, chromosome, start, end, strand, length, 
    protein ID, and gene ID, writing them to a TSV file.
    """
    # Input QC
    try:
        with open(gtf_file, 'r') as infile:
            header_lines = [line for line in infile if line.startswith('#')]
            if not header_lines:
                print("Warning: No header found. Ensure the file is a valid GTF format.")
    except FileNotFoundError:
        print(f"Error: File '{gtf_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)

    pattern = re.compile(
        r'gene_id "(.*?)";.*?transcript_id "(.*?)";.*?transcript_type "(.*?)";.*?protein_id "(.*?)"'
    )

    with open(gtf_file, 'r') as infile, open(output_tsv, 'w') as outfile:
        outfile.write("transcript_ID\tchromosome\tstart\tend\tstrand\tlength\ttranscript_type\tprotein_ID\tgene_ID\n")
        for line in infile:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'transcript':
                continue
            
            chrom, start, end, strand, attributes = fields[0], int(fields[3]), int(fields[4]), fields[6], fields[8]
            length = end - start
            match = pattern.search(attributes)
            if match:
                gene_id, transcript_id, transcript_type, protein_id = match.groups()
                protein_id = protein_id if protein_id else 'NA'
                outfile.write(f"{transcript_id}\t{chrom}\t{start}\t{end}\t{strand}\t{length}\t{transcript_type}\t{protein_id}\t{gene_id}\n")

    print(f"Parsing complete. Output saved to {output_tsv}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python parse_gtf.py <input.gtf> <output.tsv>")
        sys.exit(1)
    
    input_gtf = sys.argv[1]
    output_tsv = sys.argv[2]
    parse_gtf(input_gtf, output_tsv)
