#! bin/bash

# Take user arguments: input-, config- and output directories
while getopts i:c:o: flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        c) config_dir=${OPTARG};;
        o) output_dir=${OPTARG};;
    esac
done