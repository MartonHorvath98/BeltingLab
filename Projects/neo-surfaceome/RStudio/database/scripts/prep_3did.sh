#!/bin/bash
# Using getopt

#####################################################
## Marton Horvath, MAR 2025 - last modified: MAR 2025
#####################################################
# Usage function
usage() {
    echo "Usage: $0 -i <input_file>"
    exit 1
}

# Parse command-line options
while getopts "i:" opt; do
    case "$opt" in
        i) input_file="$(realpath $OPTARG)" ;;
        *) usage ;;
    esac
done

# Ensure input file is provided
if [ -z "$input_file" ]; then
    echo "ERROR: No input file specified."
    usage
fi

# Ensure input file exists
if [ ! -f "$input_file" ]; then
    echo "ERROR: Input file '$input_file' not found!"
    exit 1
fi

# Save input path
input_dir="$(dirname $input_file)"
# Define output names
preprocessed_file="3did_pfamInteractions.tsv"
final_output="3did_tbl.tsv"

# Save the output file paths in the same directory where the input file was found

# Step 1: Filter out PFAM domain-domain interaction information
echo "#############################################"
echo "# Step 1: Processing Pfam interactions...   #"
echo "#############################################"
less "$input_file" | grep "^#=ID" | cut -f4,5 | perl -ane '
    $F[0] =~ s/.*(PF\d+).*/$1/;
    $F[1] =~ s/.*(PF\d+).*/$1/;
    print "$F[0]\t$F[1]\n$F[1]\t$F[0]\n";
' | sort -u > "$input_dir/$preprocessed_file"

# Step 2: Process interactions into final DIMA formatted table. This is a legacy step.
echo "#############################################"
echo "# Step 2: Reformatting to DIMA format...    #"
echo "#############################################"

less "$input_dir/$preprocessed_file" | perl -ane '
    BEGIN { my %ipfam; }
    $ipfam{$F[0]}{$F[1]} = 1;
    END {
        foreach my $p1 (sort keys %ipfam) {
            foreach my $p2 (sort keys %{$ipfam{$p1}}) {
                print "$p1\t$p2\t1\t\t1.0\t\t\t\t\t\n";
            }
        }
    }
' > "$input_dir/$final_output"

echo "Processing complete! Output saved to: $input_dir/$final_output"
