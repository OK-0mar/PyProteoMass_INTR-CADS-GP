# -*- coding: utf-8 -*-
"""
Title: INTR-CADS Graduation Project
Project Title: MW of genome of E.coli
Created on Mon Aug  3 14:31:23 2026
@author: Omar Khaled
"""
#%% Import

import pandas as pd
from Bio import SeqIO
#%% Function

def calculate_mw(sequence):
    
    """
    Calculates molecular weight:
    Sum of free amino acid weights - Water_MW * (Number of residues - 1)
    """
    aa_weights = {
        'A': 89.094, 'R': 174.203, 'N': 132.119, 'D': 133.104, 'C': 121.154,
        'E': 147.130, 'Q': 146.146, 'G': 75.067,  'H': 155.156, 'I': 131.175,
        'L': 131.175, 'K': 146.189, 'M': 149.211, 'F': 165.192, 'P': 115.132,
        'S': 105.093, 'T': 119.119, 'W': 204.228, 'Y': 181.191, 'V': 117.148
    }
    water_mw = 18.015

  # Filter for standard amino acids
    valid_residues = [aa for aa in sequence.upper() if aa in aa_weights]
    num_residues = len(valid_residues)
    
    sum_of_residues = sum(aa_weights[aa] for aa in valid_residues)
    total_mw = sum_of_residues - (water_mw * (num_residues - 1))
    
    return round(total_mw, 2)

def process_genome_files(fasta_file, tabular_file, output_csv, is_ptt=False):
    
    """
    Reads FASTA via Biopython, reads PTT/CSV/Feature Table via Pandas, 
    and outputs a 2-column CSV (Locus Tag, Molecular Weight).
    """
    # 1. READ FASTA FILE USING BIOPYTHON
    print("Reading FASTA file using Biopython...")
    protein_mw_dict = {}
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_id = record.id.split()[0]
        sequence_str = str(record.seq)
        mw = calculate_mw(sequence_str)
        protein_mw_dict[protein_id] = mw
        
    print(f"Calculated MW for {len(protein_mw_dict)} proteins.")

    # 2. READ TABULAR FILE USING PANDAS
    print("\nReading tabular annotation file using Pandas...")
    
    # Try reading tab-separated first (NCBI standard), fallback to comma if needed
    try:
        df = pd.read_csv(tabular_file, sep='\t', comment='#')
        if len(df.columns) <= 1: # If it read as a single column, try comma delimiter
            df = pd.read_csv(tabular_file, sep=',')
    except Exception:
        df = pd.read_csv(tabular_file, sep=',')

    # Clean whitespace from column names (e.g. ' locus_tag ' -> 'locus_tag')
    df.columns = df.columns.str.strip()

    # Dynamic search for Protein ID column
    protein_candidates = ['Protein accession', 'product_accession', 'protein_accession', 'Protein_ID', 'PID']
    protein_col = next((col for col in protein_candidates if col in df.columns), None)

    # Dynamic search for Locus Tag column (added 'Locus tag')
    locus_candidates = ['Locus tag', 'locus_tag', 'Locus_Tag', 'Locus Tag', 'Synonym', 'locus_tag_name']
    locus_col = next((col for col in locus_candidates if col in df.columns), None)

    
    # 3. MAP MOLECULAR WEIGHTS & EXPORT RESULTS
    print("\nMapping Molecular Weights to Locus Tags...")
    df[protein_col] = df[protein_col].astype(str)
    df['Molecular_Weight_Da'] = df[protein_col].map(protein_mw_dict)
    
    # Filter down to the required two columns and drop unmapped rows
    final_df = df[[locus_col, 'Molecular_Weight_Da']].dropna()
    
    # Save to CSV using Pandas
    final_df.to_csv(output_csv, index=False)
    print(f"\nSuccess! Output saved to: {output_csv}")
#%% Execution

if __name__ == "__main__":
    
    # Input file paths
    FASTA_FILE = r"E:\INTR-CADS\Python\Ecoli_genome.faa"
    TABULAR_FILE = r"E:\INTR-CADS\Python\Ecoli_Ann.tsv"
    OUTPUT_FILE = r"E:\INTR-CADS\Python\Ecoli_locus_mw.csv"
    
    process_genome_files(FASTA_FILE, TABULAR_FILE, OUTPUT_FILE, is_ptt=False)