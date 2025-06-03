#!/usr/bin/env python3

import sys

def clean_string_ids(input_file, output_file, columns_to_clean, separator='\t'):
    columns_to_clean = set(map(int, columns_to_clean.split(',')))
    print(f"Cleaning columns: {columns_to_clean}")

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        header = next(infile) # Read the header line
        outfile.write(header)

        for line in infile:
            fields = line.strip().split(separator)

            for idx in columns_to_clean:
                if idx < len(fields) and fields[idx].startswith('9606.'):
                    fields[idx] = fields[idx].replace('9606.', '', 1)

            outfile.write('\t'.join(fields) + '\n')

    print(f"Cleaned file written to: {output_file}")


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python clean_string_ids.py <input.tsv> <output.tsv> <column_indices> <separator>")
        print("Example: python clean_string_ids.py input.tsv output.tsv 0 '\t'")
        print("         python clean_string_ids.py input.tsv output.tsv 0,2 ' '")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    column_indices = sys.argv[3]
    separator = sys.argv[4]

    clean_string_ids(input_file, output_file, column_indices, separator='\t')
