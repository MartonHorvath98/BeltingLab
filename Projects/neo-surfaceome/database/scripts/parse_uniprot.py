#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parse_uniprot.py - A script to parse the UniprotID mapping .dat file, filtering needed annotations into TSV file output.
"""

import csv
import sys
from collections import defaultdict


def parse_uniprot(input_file, output_file):
    # Structure to hold info by UniProt ID
    uniprot_data = defaultdict(lambda: {"uniprot_name": None, "preferred_name": None, "protein_ids": []})

    # Read the file and collect relevant data
    with open(input_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if len(row) < 3:
                continue  # skip malformed lines
            uniprot_id, category, value = row[0], row[1], row[2]

            if category == "UniProtKB-ID":
                uniprot_data[uniprot_id]["uniprot_name"] = value
            elif category == "Gene_Name":
                uniprot_data[uniprot_id]["preferred_name"] = value
            elif category == "STRING":
                cleaned_value = value.removeprefix("9606.")
                uniprot_data[uniprot_id]["protein_ids"].append(cleaned_value)

    # Write to output
    with open(output_file, "w", newline="") as f_out:
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(["uniprot_id", "uniprot_name", "protein_id", "preferred_name"])

        for uniprot_id, info in uniprot_data.items():
            for protein_id in info["protein_ids"]:
                writer.writerow([
                    uniprot_id,
                    info["uniprot_name"] or "",
                    protein_id,
                    info["preferred_name"] or ""
                ])

    print(f"Parsing complete. Output saved to {output_file}")

if __name__ == "__main__":    
    if len(sys.argv) != 3:
        print("Usage: python parse_uniprot.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    parse_uniprot(input_file, output_file)