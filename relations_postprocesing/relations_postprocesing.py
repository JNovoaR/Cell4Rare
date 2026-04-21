# --- INPUT DESCRIPTION --- #
# sys.argv[1]: (TSV) The "oli" file containing raw statistical results.
# sys.argv[2]: (Path) Directory containing "ngenes_" files (mapping HPOs to gene counts).
# sys.argv[3]: (String) Dataset mode: 'hpa' (Human Protein Atlas), 'ts' (Tabula Sapiens), or 'logfc' (Fold Change Differential Expression baseline relations).
# -------------------------- #

import sys
import os
import numpy as np
from statsmodels.stats.multitest import multipletests
import networkx
import obonet

# --- INITIAL CONFIGURATION --- #

f_oli = sys.argv[1]
folder_ngenes = sys.argv[2]
ts_or_hpa = sys.argv[3]

# Define HPO column index and validation file paths based on the dataset source
hpo_col = "NA"
if ts_or_hpa == "hpa":
    hpo_col = 1
    f_cell2coment = "../validation/hpa2coment"
elif ts_or_hpa == "ts":
    hpo_col = 0
    f_cell2coment = "../validation/ts2coment"
elif ts_or_hpa == "logfc":
    hpo_col = 1
    f_cell2coment = "../validation/hpa2coment"
else:
    print("Error: argv[3] must be 'hpa', 'ts', or 'logfc'", file=sys.stderr)
    exit()

# Threshold for minimum genes associated with an HPO to keep the result
maxngenes_threshold = 10
# Flag: If True, HPOs not present in the CoMent database are marked 'NA'
sensitive_2_nonexistant_in_coment = True 
f_raw_coment = "../../coment_net/DATA/ALL_CL-HPO_s"

# Load the HPO ontology graph structure
hpo_graph = obonet.read_obo("../DATA/HPO/hp.obo") 
topparent_term = "HP:0000118" # Root term: "Phenotypic abnormality"

# Define relative column indices for data processing
ngenes_col = hpo_col + 2
cluster_lvl_col = hpo_col + 3
cluster_col = hpo_col + 4
pval_col = hpo_col + 8

# --- TISSUE LIST EXTRACTION --- #
# Identifies all tissues present in the input file to locate corresponding ngenes files
l_tissues = []
with open(f_oli, "r") as f:
    for line in f:
        l_line = line.strip().split("\t")
        if ts_or_hpa in ["hpa", "logfc"]:
            tissue = l_line[hpo_col - 1]
        elif ts_or_hpa == "ts":
            tissue = l_line[cluster_col].split("/")[0]
        
        if tissue not in l_tissues:
            l_tissues.append(tissue)

# --- MAXGENES INTEGRATION --- #
# Loads gene counts per HPO/Tissue and filters the dataset by the threshold
d_hpo2ngenes = {}
for tissue in l_tissues:
    d_hpo2ngenes[tissue] = {}
    # Formatting filename based on dataset convention
    if ts_or_hpa == "ts":
        ngenes_filename = "/ngenes_" + tissue.title().replace(" ", "_")
    else:
        ngenes_filename = "/ngenes_" + tissue.lower().replace(" ", "_")
    
    with open(folder_ngenes + ngenes_filename, "r") as f:
        for line in f:
            l_line = line.strip().split("\t")
            hpo = l_line[0]
            # Take the maximum gene count reported across cells for this HPO
            max_genes = max([int(x) for x in l_line[1:]])
            d_hpo2ngenes[tissue][hpo] = max_genes

l_f_mg = [] # List to store lines meeting the gene count threshold
with open(f_oli, "r") as f:
    for line in f:
        l_line = line.strip().split("\t")
        hpo = l_line[hpo_col]
        
        if ts_or_hpa in ["hpa", "logfc"]:
            tissue = l_line[hpo_col - 1]
        elif ts_or_hpa == "ts":
            tissue = l_line[cluster_col].split("/")[0]
            
        # Check if the HPO meets the gene count threshold for this specific tissue
        if d_hpo2ngenes[tissue][hpo] >= maxngenes_threshold:
            # Reconstruct the line with the updated max_ngenes value
            l_f_mg.append("\t".join(l_
