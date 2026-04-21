# --- INPUT DESCRIPTION --- #
# This script requires multiple command-line arguments:
# sys.argv[1]: (TSV) Cell annotation file. Columns: ID, Tissue, Group, Cell_Type, Cluster.
# sys.argv[2]: (TSV) Filtered HPO file. Column 0: ID, Column 1: Name, Column 3: Gene Count.
# sys.argv[3:]: (TSVs) One or more ngenes files (from HPA2ngenes.py)
#                - Line 1: Tissue|all_genes header followed by total counts per cell.
#                - Subsequent lines: HPO_ID followed by specific counts per cell.
# -------------------------- #

import sys
import numpy as np
import math
from sklearn.linear_model import LinearRegression
from scipy.stats import kstest

# --- HELPER FUNCTIONS --- #

def convert_to_integers(string_list):
    """Converts a list of strings (from a TSV) into integers, handling errors."""
    try:
        # Standard list comprehension to cast numeric strings to actual integers
        return [int(item) for item in string_list]
    except ValueError as e:
        # Catches cases where the TSV might contain non-numeric data in data columns
        print(f"Error converting list to integers: {e}")
        return None

def cal_effect_size(Z, n1, n2):
    """
    Calculates the effect size for the Kolmogorov-Smirnov (KS) test.
    Corrects the statistic based on the two sample sizes (group vs background).
    """
    # Formula: Z / sqrt((n1*n2)/(n1+n2)) 
    # This standardizes the KS statistic to make it comparable across different group sizes
    effect_size = Z / math.sqrt((n1 * n2) / (n1 + n2))
    return effect_size

# -- CONFIGURATION -- #
# Only analyze HPO terms associated with at least 10 genes
ngenes_cutoff = 10  

# --- DATA LOADING --- #

# Paths from command line arguments
f_cell_annotation = sys.argv[1] 
f_rel_hpos = sys.argv[2]        

# Read HPO terms and filter by the ngenes_cutoff defined above
l_hpos = []
with open(f_rel_hpos, "r") as f:
    for line in f:
        l_line = line.strip().split("\t")
        hpo_id, hpo_name = l_line[0], l_line[1]
        n_related_genes = int(l_line[3]) # Index 3 contains the pre-calculated gene count
        if n_related_genes >= ngenes_cutoff:
            l_hpos.append([hpo_id, hpo_name, n_related_genes])

# Load background data: total genes expressed across all provided tissue files
l_tissues = []
l_total_expresed = [] 
# Loop through all files starting from the 3rd argument
for fname in sys.argv[3:]:
    with open(fname, "r") as f:
        line = f.readline() # Only read the first line (the header line)
        l_line = line.strip().split("\t")
        # Extract tissue name from header, e.g., "Heart|All" -> "heart"
        tissue = l_line[0].split("|")[0].lower().replace(" ", "_")
        l_tissues.append(tissue)
        # Append the integer counts of total genes expressed per cell to the global background list
        l_total_expresed += convert_to_integers(l_line[1:])

# --- CELL ANNOTATION PARSING --- #
# Groups cell indices into a hierarchical dictionary structure for targeted testing

n_cell = 0 # Tracks global cell index across all processed files
d_tissue2cell = {}
d_groupntissue2cell = {}
d_celltypengroupntissue2cell = {}
d_clusterncelltypengroupntissue2cell = {}

for included_tissue in l_tissues:
    with open(f_cell_annotation, "r") as f:
        for line in f:
            # Data cleaning: handle cases where fields might contain commas inside quotes
            if "\"" in line:
                new_line, inside_quotation = "", False
                for char in line:
                    if char == "\"":
                        inside_quotation = not inside_quotation
                        char = "" # Remove the quote mark itself
                    if char == "," and inside_quotation:
                        char = ";" # Swap comma for semicolon to preserve delimiter integrity
                    new_line += char
                line = new_line

            l_line = line.strip().split("\t")
            # Map columns to variables; uses .lower() for case-insensitive matching
            cell_id, tissue, group, cell_type, cluster = [x.lower() for x in l_line[:5]]

            # If the cell belongs to the current tissue being processed, index its metadata
            if tissue == included_tissue:
                # setdefault initializes a list if the key is new, then appends the global cell index
                d_tissue2cell.setdefault(tissue, []).append(n_cell)
                d_groupntissue2cell.setdefault(f"{tissue}/{group}", []).append(n_cell)
                d_celltypengroupntissue2cell.setdefault(f"{tissue}/{group}/{cell_type}", []).append(n_cell)
                d_clusterncelltypengroupntissue2cell.setdefault(f"{tissue}/{group}/{cell_type}/{cluster}", []).append(n_cell)
                n_cell += 1

# --- ANALYSIS LOOP --- #

n_hpos = len(l_hpos)
for n, id_n_name in enumerate(l_hpos, 1):
    hpo_id, hpo_name, n_related_genes = id_n_name
    l_related_expresed = [] # Stores HPO-
