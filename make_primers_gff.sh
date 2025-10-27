#!/bin/bash
#
# Primer Design Pipeline for GFF-based Feature Extraction
# 
# This script designs PCR primers for genomic features (e.g., exons) extracted from GFF files.
# It integrates Primer3 for primer design and MFEprimer for specificity validation.
#
# Author: Scott T Small
# Version: 1.0
# Date: October 25 2025

set -euo pipefail

show_help() {
    cat << EOF
PRIMER DESIGN PIPELINE

USAGE:
    $0 -n NAME -g GFF -r REF -f FEATURE [OPTIONS]

DESCRIPTION:
    Automated PCR primer design pipeline that:
    1. Extracts genomic features from GFF files
    2. Designs primers using Primer3
    3. Validates primer specificity using MFEprimer
    4. Generates comprehensive analysis reports

REQUIRED OPTIONS:
    -n, --name NAME         Project name/identifier for output files
    -g, --gff GFF          Path to GFF/GTF annotation file
    -r, --ref REF          Path to reference genome FASTA file  
    -f, --feature FEATURE  Feature type to extract (e.g., 'exon', 'gene', 'CDS')
    -c, --chrom CHROM      Specific chromosome to process (default: all chromosomes)
    -p, --primer3 PATH     Path to primer3_core executable 
                           (default: path_to_primer3_exe/primer3/src/primer3_core)
                           Settings file auto-detected from: {primer3_base}/settings_files/
    -m, --mfeprimer PATH   Path to mfeprimer executable
                           (default: path_to_mfeprimer_exe/mfeprimer)
    -h, --help             Show this help message and exit

EXAMPLES:

    # Process only chromosome 1 exons with custom tool paths
    $0 -n chr1_exons -g annotations.gff -r genome.fa -f exon -c chr1 -p ~/tools/primer3_core -m ~/tools/mfeprimer


REQUIREMENTS:
    Software Dependencies:
    - Primer3 (specify path with -p/--primer3)
    - MFEprimer (specify path with -m/--mfeprimer) 
    - bedtools
    - Python 3 with pandas, seaborn, matplotlib, biopython
    - GNU parallel (optional, improves performance)

    Input Files:
    - GFF/GTF annotation file with target features
    - Reference genome in FASTA format
    - primer3_v1_1_4_default_settings.txt (auto-detected from primer3 path)

OUTPUT FILES:
    - \${NAME}.\${CHROM}.\${FEATURE}.results.final.csv - Final primer results
    - primer_design.results.csv - Intermediate sorted results
    - Various diagnostic plots (PDF/PNG format)
    - mfeprimer_tmp/ - MFEprimer index and intermediate files

EOF
}

# Default values
NAME=""
GFF=""
REF=""
FEAT=""
CHROM=""
primer3_path=""
mfeprimer_path=""
SHOW_HELP=false

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--name)
            NAME="$2"
            shift 2
            ;;
        -g|--gff)
            GFF="$2"
            shift 2
            ;;
        -r|--ref)
            REF="$2"
            shift 2
            ;;
        -f|--feature)
            FEAT="$2"
            shift 2
            ;;
        -c|--chrom)
            CHROM="$2"
            shift 2
            ;;
        -p|--primer3)
            primer3_path="$2"
            shift 2
            ;;
        -m|--mfeprimer)
            mfeprimer_path="$2"
            shift 2
            ;;
        -h|--help)
            SHOW_HELP=true
            shift
            ;;
        *)
            echo "Error: Unknown option '$1'"
            echo "Use -h or --help for usage information."
            exit 1
            ;;
    esac
done

# Show help if requested
if [[ "$SHOW_HELP" == true ]]; then
    show_help
    exit 0
fi

# Validate required arguments
if [[ -z "$NAME" || -z "$GFF" || -z "$REF" || -z "$FEAT" ]]; then
    echo "Error: Missing required arguments."
    echo "Required: -n/--name, -g/--gff, -r/--ref, -f/--feature"
    echo "Use -h or --help for usage information."
    exit 1
fi

# Auto-detect primer3 settings file based on executable path
if [[ "$primer3_path" != "path_to_primer3_exe/primer3/src/primer3_core" ]]; then
    # Extract the base directory (remove /src/primer3_core)
    primer3_base_dir=$(dirname $(dirname "$primer3_path"))
    primer3_settings_file="${primer3_base_dir}/settings_files/primer3_v1_1_4_default_settings.txt"
    
    if [[ -f "$primer3_settings_file" ]]; then
        echo "Using primer3 settings file: $primer3_settings_file"
    else
        echo "Warning: Expected primer3 settings file not found at: $primer3_settings_file"
        echo "Looking for primer3_v1_1_4_default_settings.txt in current directory..."
        if [[ -f "primer3_v1_1_4_default_settings.txt" ]]; then
            primer3_settings_file="primer3_v1_1_4_default_settings.txt"
            echo "Using settings file from current directory: $primer3_settings_file"
        else
            echo "Error: primer3_v1_1_4_default_settings.txt not found!"
            echo "Please ensure the settings file exists in one of these locations:"
            echo "  1. ${primer3_base_dir}/settings_files/"
            echo "  2. Current working directory"
            exit 1
        fi
    fi
else
    # Using default path - check current directory
    primer3_settings_file="primer3_v1_1_4_default_settings.txt"
    if [[ ! -f "$primer3_settings_file" ]]; then
        echo "Error: primer3_v1_1_4_default_settings.txt not found in current directory"
        echo "Please provide the path to primer3 executable with -p/--primer3 flag"
        echo "or place the settings file in the current directory"
        exit 1
    fi
fi

if [[ ! -f $GFF ]]; then
    echo "Error: GFF file $GFF not found."
    exit 1
fi

if [[ ! -f $REF ]]; then
    echo "Error: Reference FASTA $REF not found."
    exit 1
fi

# Setup MFEprimer index once (check if already exists)
IDX_DIR="mfeprimer_tmp/${NAME}_idx"
REF_BASENAME=$(basename "$REF")

if [[ ! -d "$IDX_DIR" ]] || [[ ! -f "$IDX_DIR/${REF_BASENAME}.primerqc" ]]; then
    echo "Setting up MFEprimer index for ${NAME}..."
    mkdir -p "$IDX_DIR"
    # Create absolute path for reference
    REF_ABS=$(realpath "$REF")
    ln -sf "$REF_ABS" "$IDX_DIR/"
    (
        cd "$IDX_DIR"
        $mfeprimer_path index -i "${REF_BASENAME}"
    )
    echo "MFEprimer index created."
else
    echo "MFEprimer index already exists for ${NAME}, skipping creation."
fi

# Get list of chromosomes if not specified
if [[ -z "$CHROM" ]]; then
    mapfile -t CHROMS < <(cut -f1 $GFF | sort -u)
else
    CHROMS=($CHROM)
fi

for CHR in "${CHROMS[@]}"; do
    echo "Processing chromosome: $CHR"
    grep -P "^$CHR\t.*\t$FEAT\t" $GFF | cut -f1,4-5,9 > ${NAME}.${FEAT}.${CHR}.txt || continue
    awk '{ print $1"\t"$2-1"\t"$3"\t"$4; }' ${NAME}.${FEAT}.${CHR}.txt > ${NAME}.${FEAT}.${CHR}.bed
    bedtools getfasta -fi $REF -fo ${NAME}.${FEAT}.${CHR}.fasta -bed ${NAME}.${FEAT}.${CHR}.bed

    # Run primer3
    python primer_pipeline.py primer3-input ${NAME}.${FEAT}.${CHR}.fasta ${NAME}.${FEAT}.${CHR}.primer3_input.dat
    $primer3_path ${NAME}.${FEAT}.${CHR}.primer3_input.dat --p3_settings_file="$primer3_settings_file" > ${NAME}.${FEAT}.${CHR}.primer3_output.dat

    # Run MFEprimer
    FASTA_DIR="mfeprimer_tmp/fasta_in/${CHR}"
    mkdir -p "$FASTA_DIR"
    
    # Copy primer3 output to fasta directory
    cp "${NAME}.${FEAT}.${CHR}.primer3_output.dat" "$FASTA_DIR/"
    
    # Run mfeprimer_input.py in the fasta directory
    (
        cd "$FASTA_DIR"
        python ../../../primer_pipeline.py mfeprimer-input "${NAME}.${FEAT}.${CHR}.primer3_output.dat" "${NAME}" || { echo "mfeprimer-input failed for $CHR"; exit 1; }
        
        # Run MFEprimer on each fasta file using absolute path to index
        IDX_ABS=$(realpath "../../${NAME}_idx")
        if command -v parallel >/dev/null 2>&1; then
            parallel "$mfeprimer_path -d ${IDX_ABS}/${REF_BASENAME} -i {} -j -o {.}.results" ::: *.fa
        else
            # Fallback if parallel is not available
            for fa_file in *.fa; do
                if [[ -f "$fa_file" ]]; then
                    output_name="${fa_file%.*}.results"
                    $mfeprimer_path -d "${IDX_ABS}/${REF_BASENAME}" -i "$fa_file" -j -o "$output_name"
                fi
            done
        fi
    )
done

# Now run python scripts to sort and analyze primers
echo "Sorting and analyzing primers..."
if [[ -z "$CHROM" ]]; then
    CHROM_STR="all_chroms"
else
    CHROM_STR="$CHROM"
fi

# Sort and analyze results
python primer_pipeline.py sort-analyze --output "${NAME}.${CHROM_STR}.${FEAT}.results.final.csv"

# Verify primers
echo "Verifying primer results..."
python primer_pipeline.py verify

echo "Primer design pipeline completed successfully!"
echo "Final results: ${NAME}.${CHROM_STR}.${FEAT}.results.final.csv"
