# --- INPUT DESCRIPTION --- #
# sys.argv[1]: (TSV) HPO-to-Gene annotation file. (phenotype_to_genes.txt directly downloaded from HPO page)
#               Expected columns: [HPO_ID, HPO_Name, Entrez_ID, Gene_Symbol, ...]
# sys.argv[2]: (TSV) ID Mapping file. 
#               Expected columns: [Ensembl_ID, Entrez_ID]
# Output:       Saves 'hpo2genes.json' in the current working directory.
# -------------------------- #

import json
from datetime import datetime
import sys

# --- DATA LOADING --- #

# Path to the file containing Ensembl-to-Entrez ID relationships
ensembl_entrez_map = sys.argv[2]

# Dictionary to store Entrez IDs as keys and Ensembl IDs as values for fast lookup
entrez_ensembl = {}
with open(ensembl_entrez_map) as f: 
    for line in f:
        # Standard TSV parsing: remove newline and split by tab
        l_line = line.strip("\n").split("\t")
        # Map Entrez (column 1) to Ensembl (column 0)
        entrez_ensembl[l_line[1]] = l_line[0]

# Load the HPO phenotype-to-gene annotation file
# Read all lines, skip the header ([1:]), and store as a list of lists
with open(sys.argv[1], 'r') as f: 
    phe_to_genes = [line.strip('\n').split('\t') for line in f.readlines()][1:]


# --- HPO-GENE RELATIONSHIP MAPPING --- #

# Initialize a dictionary where keys are HPO IDs and values are metadata structures
hpo_genes = {
    line[0]: {'HPOname': line[1], 'EntrezID': [], 'EnsemblID': [], 'GeneName': []} 
    for line in phe_to_genes
}

# Fill the metadata structures by checking if an Ensembl translation exists
for line in phe_to_genes:
    hpo_id = line[0]
    entrez_id = line[2]
    gene_name = line[3]
    
    # Only process rows where the Entrez ID has a corresponding Ensembl ID
    if entrez_id in entrez_ensembl.keys():
        hpo_genes[hpo_id]['EntrezID'].append(entrez_id)
        # Perform the ID translation
        hpo_genes[hpo_id]['EnsemblID'].append(entrez_ensembl[entrez_id])
        hpo_genes[hpo_id]['GeneName'].append(gene_name)


# --- CLEANUP (OPTIONAL) --- #
# The following block (commented out) would remove HPO terms that have zero mapped genes.
## hp_to_remove = []
## for hp in hpo_genes:
##     if len(hpo_genes[hp]['EntrezID']) == 0: 
##         hp_to_remove.append(hp)
##
## for hp in hp_to_remove: 
##     del hpo_genes[hp]


# --- EXPORT DATA --- #

# Write the final dictionary to a JSON file with 4-space indentation for readability
with open('hpo2genes.json', 'w') as f: 
    json.dump(hpo_genes, f, indent=4)
