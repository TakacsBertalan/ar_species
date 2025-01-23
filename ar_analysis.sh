#!/bin/bash

# Check if -d flag is provided
if [[ "$1" == "-d" && -n "$2" ]]; then
    WORKDIR="$2"
    shift 2
else
    echo "Please provide an input folder"
    echo "Usage: run_microbiome -d path/to/folder"
    exit 1
fi

snakemake -d "$WORKDIR" -s /media/deltagene/microbiome_2/AR_snakemake/Snakefile "$@"
