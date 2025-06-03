#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re

def parse_hmmscan_output(input_file, assoc_output, domain_output):
    domain_seen = {}
    
    with open(input_file, 'r') as infile, \
         open(assoc_output, 'w') as assoc_out, \
         open(domain_output, 'w') as domain_out:
        
        # Headers
        assoc_out.write("transcript_id\tprotein_length\tdomain_id\tdomain_length\tdomain_start\tdomain_end\taccuracy\tscore\tE-value\n")
        domain_out.write("domain_id\tdomain_name\tdescription\n")


        for line in infile:
            if line.startswith('#') or line.strip() == '':
                continue
            
            # Split line into tokens, preserve the final field (description)
            parts = re.split(r'\s{1,}', line.strip(), maxsplit=22)
            if len(parts) < 23:
                continue  # skip malformed lines

            domain_name = parts[0]                 # "target name"
            domain_id = parts[1].split('.')[0]     # "accession"
            domain_length = parts[2]               # "tlen"
            transcript_id = parts[3].split('.')[0] # "query name"
            transcript_length = parts[5]           # "qlen"
            e_value = parts[12]                     # "E-value"
            domain_score = parts[13]               # "(this domain) score"
            hmm_start = parts[15]                  # "(hmm coord) from"
            hmm_end = parts[16]                    # "(hmm coord) to"
            acc = parts[21]                        # "acc"
            description = parts[22]                # "description of target"

             # Write to domain table if not seen
            if domain_id not in domain_seen:
                domain_out.write(f"{domain_id}\t{domain_name}\t{description}\n")
                domain_seen[domain_id] = True

            # Write to association table
            assoc_out.write(f"{transcript_id}\t{transcript_length}\t{domain_id}\t{domain_length}\t{hmm_start}\t{hmm_end}\t{acc}\t{domain_score}\t{e_value}\n")

    print(f"Parsing complete.\n- Transcript-domain associations saved to '{assoc_output}'\n- Unique domain records saved to '{domain_output}'.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python parse_domains.py <input_file> <assoc_output.tsv> <domain_output.tsv>")
        sys.exit(1)

    parse_hmmscan_output(sys.argv[1], sys.argv[2], sys.argv[3])
