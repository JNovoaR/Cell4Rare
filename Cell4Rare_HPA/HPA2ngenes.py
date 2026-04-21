# --- INPUT DESCRIPTION --- #
# This script requires 4 command-line arguments:
# sys.argv[1]: (TSV) Gene expression matrix. Rows = Cells/Samples, Cols = Genes.
# sys.argv[2]: (TSV) List of HPO terms to analyze. Column 0: ID, Column 1: Name.
# sys.argv[3]: (JSON) HPO to genes mapping file where keys are HPO IDs and values contain 'EnsemblID' lists.
# sys.argv[4]: (String) Name of the tissue being processed (used for output headers).
# -------------------------- #

import sys
import numpy as np
import json
import math
import time
from datetime import datetime
from scipy.sparse import csr_matrix

#-- FIXED PARAMS --#
# Minimum number of genes an HPO term must be associated with to be included
ngenes_cutoff = 10 
# Toggle to print performance metrics to stderr
timing = True
#------------------#

def tsv_to_csr_dense(tsv_filepath):
    """
    Reads a gene expression TSV, skips headers/labels, and converts 
    it into a Compressed Sparse Row (CSR) matrix for memory efficiency.
    """
    try:
        # Load numeric data, skipping the first row (header) and first column (gene/cell names)
        # Uses [:,1:] slicing to discard the index column
        dense_matrix = np.loadtxt(tsv_filepath, delimiter='\t', skiprows=1, dtype=float)[:,1:]
        # Convert the dense numpy array into a sparse CSR matrix to save memory on zero-values
        csr_matrix_result = csr_matrix(dense_matrix)
        return csr_matrix_result
    except Exception as e:
        # Catch errors like FileNotFoundError or parsing issues
        print(f"Error reading or converting TSV: {e}")
        return None

def sparse_to_dense(sparse_list, default_value=0):
    """
    Expands a list of (index, value) tuples into a full-length dense list.
    """
    if not sparse_list:
        return []
    # Identify the required list length based on the highest index provided
    max_index = sparse_list[-1][0]  
    # Pre-fill a list with the default value (usually 0)
    dense_list = [default_value] * (max_index + 1)

    # Map the sparse values into their specific positions in the dense list
    for index, value in sparse_list:
        dense_list[index] = value

    return dense_list

# --- INITIALIZATION --- #

# Log the script start time to stderr for tracking execution logs
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print("Launched at: =", dt_string, file=sys.stderr)

# Load the HPA (Human Protein Atlas) expression file from the first command line argument
HPA_file = sys.argv[1]
raw_counts = tsv_to_csr_dense(HPA_file)

# Extract gene names from the header of the HPA file to create a mapping index
with open(HPA_file, "r") as f:
    for line in f:
        l_line = line.strip().split("\t")
        # Store all columns except the first one (which identifies the row/sample)
        l_genes_wvariants = l_line[1:] 
        break # Only read the first line (header)

# Clean gene IDs: remove version numbers (e.g., "ENSG0000.1" -> "ENSG0000") for easier matching
l_genes = [gene.split(".")[0] for gene in l_genes_wvariants]

# Get the tissue name from the fourth command line argument
tissue_name = sys.argv[4].lower()

# Load the HPO-to-Gene mapping JSON provided in the third argument
with open(sys.argv[3], 'r') as f:
    hpo2genes = json.load(f)

# Load the list of HPO terms to study from the second argument
l_hpos = []
with open(sys.argv[2], "r") as f:
    for line in f:
        l_line = line.strip().split("\t")
        hpo_id = l_line[0] # The ID (e.g., HP:0001234)
        hpo_name = l_line[1] # The human-readable name
        
        try:
            # Only keep HPO terms that have enough associated genes in the database to be statistically relevant
            n_related_genes = len(hpo2genes[hpo_id]['EnsemblID'])
            if n_related_genes >= ngenes_cutoff:
                l_hpos.append([hpo_id, hpo_name])
        except KeyError:
            # Skip HPO IDs that exist in the study list but are missing from the reference mapping
            pass 

print("INPUT FILES READ SUCCESSFULLY!", file=sys.stderr)
print("PREPROCESSING...", file=sys.stderr)

# --- STATISTICAL HELPER --- #

def cal_effect_size(Z, n1, n2):
    """
    Calculates the effect size for the Kolmogorov-Smirnov statistic,
    adjusting for the sample sizes of the two groups being compared.
    """
    # Standard formula for normalizing the statistic by sample size
    effect_size = Z * math.sqrt((n1 * n2) / (n1 + n2))
    return effect_size

# --- MAIN PROCESSING LOOP --- #

n_hpos = len(l_hpos)
print("DONE!", f"Number of HPOs to be processed: {n_hpos}", file=sys.stderr)

l_total_expresed = [] # Stores total count of expressed genes (any gene) per cell
first_lap = True      # Flag to ensure background 'total' counts are calculated only once

n = 0
# Iterate through each filtered HPO term
for hpoid_and_name in l_hpos:
    start_time = time.perf_counter() # Start timer for this specific HPO term
    hpo = hpoid_and_name[0]
    hponame = hpoid_and_name[1]
    l_related_expresed = [] # Stores count of HPO-specific genes expressed in each cell
    n += 1
    
    print(f"Processing {hpo} ({hponame}) [{n}/{n_hpos}]", file=sys.stderr)
    
    # Use a set for the current HPO's genes to allow O(1) lookup speed
    related_genes = set(hpo2genes[hpo]['EnsemblID']) 

    # Iterate through every row (cell) in the sparse matrix
    cell_indices = range(0, raw_counts.shape[0])
    for n_cell in cell_indices:
        n_related = 0
        
        # Access the indices of non-zero entries for this cell directly from the CSR structure
        non0_indices = raw_counts[n_cell].indices
        
        # Loop through indices of genes that have expression > 0
        for expressed_gene_index in non0_indices:
            expressed_gene = l_genes[expressed_gene_index]
            # Check if this expressed gene is one of the ones linked to the current HPO
            if expressed_gene in related_genes:
                n_related += 1
        
        # On the very first HPO loop, build the background reference of total active genes per cell
        if first_lap:
            l_total_expresed.append(len(non0_indices))
        
        # Track how many HPO-relevant genes were found in this cell
        l_related_expresed.append(n_related)

    # --- OUTPUT RESULTS --- #
    
    # Print the background total counts header once at the start of the output stream
    if first_lap:
        print(f"{tissue_name}|all_genes\t" + "\t".join(str(x) for x in l_total_expresed))
    
    # Print the specific HPO ID followed by the tab-separated counts for every cell
    print(f"{hpo}\t" + "\t".join(str(x) for x in l_related_expresed))
    
    # Flip flag so we don't recalculate total counts or reprint header in next iteration
    first_lap = False 
    
    # Output performance time for the current HPO if timing is enabled
    if timing:
        end_time = time.perf_counter()
        print(f"Elapsed time: {end_time - start_time:.4f} seconds", file=sys.stderr)
