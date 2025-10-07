#!/usr/bin/env python3
"""
Primer Design Pipeline - Consolidated Python Module

This module contains all the Python functionality for the primer design pipeline:
- Primer3 input preparation
- MFEprimer input preparation  
- Results sorting and analysis
- Primer verification

Author: Consolidated primer design pipeline
Version: 2.0
"""

import json
from pathlib import Path
import argparse
import sys
from glob import glob

# Optional imports with graceful fallbacks
try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False
    print("Warning: pandas not available - some functionality will be limited")

try:
    import seaborn as sns
    import matplotlib.pyplot as plt
    HAS_PLOTTING = True
except ImportError:
    HAS_PLOTTING = False
    print("Warning: seaborn/matplotlib not available - plots will be skipped")

try:
    from Bio import SeqIO
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("Warning: BioPython not available - primer3-input command not available")


def prepare_primer3_input(fasta_file, output_file):
    """
    Convert FASTA sequences to Primer3 input format.
    
    Args:
        fasta_file (str): Input FASTA file path
        output_file (str): Output Primer3 input file path
    """
    if not HAS_BIOPYTHON:
        print("Error: BioPython is required for primer3-input command")
        sys.exit(1)
        
    print(f"Preparing Primer3 input: {fasta_file} -> {output_file}")
    
    with open(output_file, 'w') as f:
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            f.write(f"SEQUENCE_ID={seq_record.id}\n")
            f.write(f"SEQUENCE_TEMPLATE={seq_record.seq}\n=\n")
    
    print(f"Primer3 input file created: {output_file}")


def prepare_mfeprimer_input(primer3_output, output_prefix):
    """
    Convert Primer3 output to MFEprimer input FASTA files.
    
    Args:
        primer3_output (str): Primer3 output file path
        output_prefix (str): Prefix for output FASTA files
    """
    if not HAS_PANDAS:
        print("Error: pandas is required for mfeprimer-input command")
        sys.exit(1)
        
    print(f"Preparing MFEprimer input from: {primer3_output}")
    
    # Parse Primer3 output
    primer_records = []
    with open(primer3_output) as f:
        record = {}
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line == "=":
                if record:
                    primer_records.append(record)
                    record = {}
                continue
            if "=" in line:
                key, value = line.split("=", 1)
                record[key] = value

    # Create DataFrame
    primer_df = pd.DataFrame({
        "search_string": [r["SEQUENCE_ID"] for r in primer_records],
        "primer3_forward": [r.get("PRIMER_LEFT_0_SEQUENCE", "") for r in primer_records],
        "primer3_reverse": [r.get("PRIMER_RIGHT_0_SEQUENCE", "") for r in primer_records]
    })

    # Filter valid primers
    primer_df = primer_df.dropna(subset=["primer3_forward", "primer3_reverse"])
    primer_df = primer_df[(primer_df["primer3_forward"] != "") & (primer_df["primer3_reverse"] != "")]

    # Create FASTA files for MFEprimer
    files_created = 0
    for _, row in primer_df.iterrows():
        search_string = row["search_string"]
        primer3_forward = row["primer3_forward"]
        primer3_reverse = row["primer3_reverse"]

        output_file = f"{output_prefix}.{search_string}.fa"
        with open(output_file, "w") as f:
            f.write(f">{search_string}_F\n{primer3_forward}\n")
            f.write(f">{search_string}_R\n{primer3_reverse}\n")
        files_created += 1
    
    print(f"Created {files_created} MFEprimer input files")


def parse_primer_json(primer_json):
    """
    Parse MFEprimer JSON output into a pandas DataFrame.
    
    Args:
        primer_json (str): Path to MFEprimer JSON result file
        
    Returns:
        pd.DataFrame: Parsed primer information
    """
    primer_json = Path(primer_json)
    name = primer_json.stem  # this replaces split and replace

    # Read in json file
    with open(primer_json) as f:
        primer_dict = json.load(f)

    primer_info_dict = {
        "file_name": [name.strip(".results")],
        "chromosome": [],
        "start": [],
        "end": [],
        "search_string_0ix": [],
        "id_forward": [],
        "primer3_forward": [],
        "id_reverse": [],
        "primer3_reverse": [],
        "Ta_amp": [],
        "Tm_diff": [],
        "amplicon_length": [],
        "amplicon_sequence": [],
        "gc_content_amplicon": [],
        "dimer_count": [],
        "hairpin_count": [],
        "snp_count_forward": [],
        "snp_count_reverse": [],
        "amplicon_count": [],
        "Ta_ampN": [],
        "length_ampN": [],
        "orient_ampN": [],
    }

    def len_or_zero(x):
        return 0 if x is None else len(x)

    def calc_gc(seq):
        gc = sum(1 for base in seq.upper() if base in "GC")
        return gc / len(seq) * 100 if len(seq) > 0 else 0

    def revcomp(seq):
        return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]

    amplist = primer_dict.get("AmpList", [])
    if amplist:
        primer_name = amplist[0].get("F").get("Seq").get("ID").replace("_F", "").replace(":", "-")
        chrom, start, stop = primer_name.split("-")
        primer_info_dict["chromosome"] = [chrom]
        primer_info_dict["start"] = [int(start.split("_")[0])]
        primer_info_dict["end"] = [int(stop.split("_")[0])]
        primer_info_dict.get("search_string_0ix").append(f"{chrom}:{int(start)}-{stop}")
        primer_info_dict.get("id_forward").append(amplist[0].get("F").get("Seq").get("ID"))
        primer_info_dict.get("primer3_forward").append(amplist[0].get("F").get("Seq").get("Seq"))
        primer_info_dict.get("id_reverse").append(amplist[0].get("R").get("Seq").get("ID"))
        primer_info_dict.get("primer3_reverse").append(revcomp(amplist[0].get("R").get("Seq").get("Seq")))
        primer_info_dict.get("Ta_amp").append(round(amplist[0].get("Ta")))
        primer_info_dict["Tm_diff"].append(abs(round(amplist[0].get("F").get("Tm") - amplist[0].get("R").get("Tm"))))
        primer_info_dict.get("amplicon_length").append(amplist[0].get("Size"))
        primer_info_dict["amplicon_sequence"].append(amplist[0].get("P").get("Seq")["Seq"])
        primer_info_dict["gc_content_amplicon"].append(round(amplist[0].get("GC")))
        primer_info_dict.get("dimer_count").append(len_or_zero(primer_dict.get("DimerList")))
        primer_info_dict.get("hairpin_count").append(len_or_zero(primer_dict.get("HairpinList")))
        primer_info_dict.get("snp_count_forward").append(amplist[0].get("F").get("SnpCount"))
        primer_info_dict.get("snp_count_reverse").append(amplist[0].get("R").get("SnpCount"))
        primer_info_dict.get("amplicon_count").append(len_or_zero(primer_dict.get("AmpList")))
        
        if primer_info_dict.get("amplicon_count")[0] > 1:
            ta_ampn = []
            len_ampn = []
            orient = []
            for i in range(1, primer_info_dict.get("amplicon_count")[0]):
                ta = amplist[i].get("Ta")
                if abs(ta - primer_info_dict["Ta_amp"][0]) < 10:
                    ta_ampn.append(round(ta))
                    len_ampn.append(amplist[i].get("Size"))
                    ampn_F = amplist[1].get("F")["Seq"]["ID"].split("_")[-1]
                    ampn_R = amplist[1].get("R")["Seq"]["ID"].split("_")[-1]
                    orient.append(ampn_F + ampn_R)
                else:
                    continue
        else:
            ta_ampn = [None]
            len_ampn = [None]
            orient = [None]
        
        primer_info_dict.get("Ta_ampN").append(ta_ampn)
        primer_info_dict.get("length_ampN").append(len_ampn)
        primer_info_dict.get("orient_ampN").append(orient)
    
    return pd.DataFrame(primer_info_dict)


def sort_and_analyze_primers(output_file="primer_design.results.final.csv"):
    """
    Sort MFEprimer results and generate analysis plots.
    
    Args:
        output_file (str): Final output file name
    """
    if not HAS_PANDAS:
        print("Error: pandas is required for sort-analyze command")
        sys.exit(1)
        
    print("Sorting and analyzing primer results...")
    
    # Path to MFEprimer output
    base_dir = Path("mfeprimer_tmp/fasta_in")
    primer_df = pd.DataFrame()
    
    if not base_dir.exists():
        print(f"Warning: MFEprimer output directory {base_dir} not found")
        return
    
    # Collect all JSON results
    json_files_found = 0
    for chrom_dir in base_dir.iterdir():
        if not chrom_dir.is_dir():
            continue
        for json_file in chrom_dir.glob("*.results.json"):
            try:
                result_df = parse_primer_json(json_file)
                primer_df = pd.concat([primer_df, result_df], ignore_index=True)
                json_files_found += 1
            except Exception as e:
                print(f"Warning: Failed to parse {json_file}: {e}")
    
    print(f"Processed {json_files_found} JSON result files")
    
    if primer_df.empty:
        print("No primer results found to analyze")
        return
    
    # Calculate additional metrics
    primer_df["masked_forward"] = primer_df["primer3_forward"].apply(lambda s: sum(1 for c in s if c.islower()))
    primer_df["masked_reverse"] = primer_df["primer3_reverse"].apply(lambda s: sum(1 for c in s if c.islower()))

    # Amplicon count distribution plot
    if HAS_PLOTTING:
        amplicon_counts = primer_df["amplicon_count"].value_counts().sort_index()
        plt.figure(figsize=(10, 5))
        sns.barplot(x=amplicon_counts.index, y=amplicon_counts.values)
        plt.xlabel("Amplicon Count")
        plt.ylabel("Frequency")
        plt.title("Distribution of Amplicon Counts")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig("amplicon_count_distribution.pdf")
        plt.close()

    # Calculate primer scoring
    primer_df.loc[:, "score"] = (
        primer_df["Tm_diff"] * 1.5 +
        abs(primer_df["gc_content_amplicon"] - 50) * 1.0 +
        primer_df["amplicon_count"] * 2 +
        primer_df["dimer_count"] * 2 +
        primer_df["hairpin_count"] * 2 + 
        primer_df["masked_forward"] / primer_df["amplicon_length"] * 100 + 
        primer_df["masked_reverse"] / primer_df["amplicon_length"] * 100
    )

    # Sort by score (lower is better)
    sorted_df = primer_df.sort_values("score")
    sorted_df.to_csv("primer_design.results.csv", index=False)

    # Diagnostic figures for single-amplicon primers only
    df = primer_df[primer_df["amplicon_count"] == 1].copy()
    print(f"Saved sorted primer list to primer_design.results.csv, found {len(df)} primers w/ no off target")
    
    if not df.empty:
        generate_diagnostic_plots(df)
    
    # Save final results
    sorted_df.to_csv(output_file, index=False)
    print(f"Final results saved to {output_file}")


def generate_diagnostic_plots(df):
    """Generate diagnostic plots for primer analysis."""
    if not HAS_PLOTTING:
        print("Skipping diagnostic plots - matplotlib/seaborn not available")
        return
    
    # Tm difference
    plt.figure()
    sns.histplot(df["Tm_diff"], bins=20, kde=True)
    plt.title("Distribution of Tm Difference (F - R)")
    plt.xlabel("Tm Difference (°C)")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig("tm_difference_distribution.png")
    plt.close()

    # Ta  
    plt.figure()
    sns.histplot(df["Ta_amp"], bins=20, kde=True)
    plt.title("Distribution of Ta")
    plt.xlabel("Annealing Temp (°C)")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig("Ta_distribution.png")
    plt.close()

    # Amplicon length
    plt.figure()
    sns.histplot(df["amplicon_length"], bins=20, kde=True)
    plt.title("Amplicon Length Distribution")
    plt.xlabel("Amplicon Length (bp)")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig("amplicon_length_distribution.png")
    plt.close()

    # Dimer count
    plt.figure()
    sns.countplot(x="dimer_count", data=df)
    plt.title("Dimer Count per Primer Pair")
    plt.xlabel("Dimer Count")
    plt.ylabel("Number of Primer Pairs")
    plt.tight_layout()
    plt.savefig("dimer_count_distribution.png")
    plt.close()

    # Hairpin count
    plt.figure()
    sns.countplot(x="hairpin_count", data=df)
    plt.title("Hairpin Count per Primer Pair")
    plt.xlabel("Hairpin Count")
    plt.ylabel("Number of Primer Pairs")
    plt.tight_layout()
    plt.savefig("hairpin_count_distribution.png")
    plt.close()

    # GC content plots
    plt.figure()
    sns.histplot(df[["gc_content_amplicon"]].melt(value_name="GC%", var_name="Amplicon"), x="GC%", bins=20, kde=True)
    plt.title("GC Content Distribution of Amplicons")
    plt.xlabel("GC Content (%)")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig("gc_content_distribution.png")
    plt.close()


def verify_primers():
    """
    Verify primer results by comparing with original Primer3 output.
    """
    if not HAS_PANDAS:
        print("Error: pandas is required for verify command")
        sys.exit(1)
        
    print("Verifying primer results...")
    
    # Parse PRIMER3 output in key=value format
    primer3_files = glob("*.primer3_output.dat")
    if not primer3_files:
        print("Warning: No .primer3_output.dat files found for verification.")
        return
    
    all_primer_records = []
    for file in primer3_files:
        primer_records = []
        with open(file) as f:
            record = {}
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line == "=":
                    if record:
                        primer_records.append(record)
                        record = {}
                    continue
                if "=" in line:
                    key, value = line.split("=", 1)
                    record[key] = value
        all_primer_records.extend(primer_records)

    if not all_primer_records:
        print("No primer records found in Primer3 output files")
        return

    primer_df = pd.DataFrame({
        "search_string_0ix": [r["SEQUENCE_ID"] for r in all_primer_records],
        "primer3_forward": [r.get("PRIMER_LEFT_0_SEQUENCE", "").upper() for r in all_primer_records],
        "primer3_reverse": [r.get("PRIMER_RIGHT_0_SEQUENCE", "").upper() for r in all_primer_records]
    })

    primer_df = primer_df.dropna(subset=["primer3_forward", "primer3_reverse"])
    primer_df = primer_df[(primer_df["primer3_forward"] != "") & (primer_df["primer3_reverse"] != "")]
    primer_df = primer_df.drop_duplicates(subset="search_string_0ix", keep="first")

    # Load annotated CSV
    if not Path("primer_design.results.csv").exists():
        print("Warning: primer_design.results.csv not found for verification")
        return
        
    annotated_df = pd.read_csv("primer_design.results.csv")

    # Merge and compare
    merged = annotated_df.merge(primer_df, on="search_string_0ix", how="left")
    merged = merged.dropna()
    
    if not merged.empty:
        merged["primer3_forward_match"] = merged["primer3_forward_x"].str.upper() == merged["primer3_forward_y"]
        merged["primer3_reverse_match"] = merged["primer3_reverse_x"].str.upper() == merged["primer3_reverse_y"]
        
        # Save results
        merged.to_csv("primer_design.results.final.csv", index=False)
        print("Verification completed. Results saved to primer_design.results.final.csv")
        
        # Report verification statistics
        forward_matches = merged["primer3_forward_match"].sum()
        reverse_matches = merged["primer3_reverse_match"].sum()
        total = len(merged)
        
        print("Verification Statistics:")
        print(f"  Total primers verified: {total}")
        print(f"  Forward primer matches: {forward_matches}/{total} ({forward_matches/total*100:.1f}%)")
        print(f"  Reverse primer matches: {reverse_matches}/{total} ({reverse_matches/total*100:.1f}%)")
    else:
        print("No primer data available for verification")


def main():
    """Main command-line interface for the primer pipeline."""
    
    parser = argparse.ArgumentParser(
        description="Primer Design Pipeline - Python Components",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Subcommands:
  primer3-input    Convert FASTA to Primer3 input format
  mfeprimer-input  Convert Primer3 output to MFEprimer input
  sort-analyze     Sort and analyze MFEprimer results  
  verify           Verify primer results against original Primer3 output
  
Examples:
  python primer_pipeline.py primer3-input sequences.fasta primer3.input
  python primer_pipeline.py mfeprimer-input primer3.output project_name
  python primer_pipeline.py sort-analyze --output final_results.csv
  python primer_pipeline.py verify
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Primer3 input subcommand
    p3_parser = subparsers.add_parser('primer3-input', help='Prepare Primer3 input from FASTA')
    p3_parser.add_argument('fasta_file', help='Input FASTA file')
    p3_parser.add_argument('output_file', help='Output Primer3 input file')
    
    # MFEprimer input subcommand  
    mfe_parser = subparsers.add_parser('mfeprimer-input', help='Prepare MFEprimer input from Primer3 output')
    mfe_parser.add_argument('primer3_output', help='Primer3 output file')
    mfe_parser.add_argument('output_prefix', help='Output prefix for FASTA files')
    
    # Sort and analyze subcommand
    sort_parser = subparsers.add_parser('sort-analyze', help='Sort and analyze primer results')
    sort_parser.add_argument('--output', default='primer_design.results.final.csv', help='Output file name')
    
    # Verify subcommand
    subparsers.add_parser('verify', help='Verify primer results')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return
    
    try:
        if args.command == 'primer3-input':
            prepare_primer3_input(args.fasta_file, args.output_file)
        elif args.command == 'mfeprimer-input':
            prepare_mfeprimer_input(args.primer3_output, args.output_prefix)
        elif args.command == 'sort-analyze':
            sort_and_analyze_primers(args.output)
        elif args.command == 'verify':
            verify_primers()
        else:
            parser.print_help()
    
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()