#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Parse GTF file and extract relevant information for the database annotation of each transcript.
"""

import sys
import re


def parse_gtf(gtf_file, output_tsv):
    """
    Parses the Ensembl GTF file and extracts transcript ID, chromosome, start, end, strand, length,
    transcript type, protein ID, gene ID, gene symbol, and HGNC symbol.
    """
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

    with open(gtf_file, 'r') as infile, open(output_tsv, 'w') as outfile:
        outfile.write("#transcript_id\tchromosome\tstart\tend\tstrand\tlength\ttranscript_type\tprotein_id\tgene_id\tgene_symbol\thgnc_symbol\n")
        for line in infile:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'transcript':
                continue

            chrom, start, end, strand, attributes = fields[0], int(fields[3]), int(fields[4]), fields[6], fields[8]
            length = end - start

            # Dictionary-based attribute extraction
            attr_dict = {}
            for attr in attributes.strip().split(';'):
                if attr.strip() == '':
                    continue
                key, val = attr.strip().split(' ', 1)
                attr_dict[key] = val.strip('"')

            # Required fields (some may be missing)
            gene_id = attr_dict.get("gene_id", "NA").split('.')[0]
            transcript_id = attr_dict.get("transcript_id", "NA").split('.')[0]
            transcript_type = attr_dict.get("transcript_type", "NA")
            protein_id = attr_dict.get("protein_id", "NA").split('.')[0] if "protein_id" in attr_dict else "NA"
            gene_symbol = attr_dict.get("gene_name", "NA")
            hgnc_symbol = attr_dict.get("hgnc_id", "NA").replace("HGNC:", "") if "hgnc_id" in attr_dict else "NA"

            outfile.write(f"{transcript_id}\t{chrom}\t{start}\t{end}\t{strand}\t{length}\t{transcript_type}\t{protein_id}\t{gene_id}\t{gene_symbol}\t{hgnc_symbol}\n")

    print(f"Parsing complete. Output saved to {output_tsv}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python parse_gtf.py <input.gtf> <output.tsv>")
        sys.exit(1)

    input_gtf = sys.argv[1]
    output_tsv = sys.argv[2]
    parse_gtf(input_gtf, output_tsv)
