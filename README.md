# ceph_primer_pipeline
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
    -n, --name NAME        Project name/identifier for output files
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
