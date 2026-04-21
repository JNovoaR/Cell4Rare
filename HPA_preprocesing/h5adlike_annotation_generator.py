# --- INPUT DESCRIPTION --- #
# sys.argv[1]: (TSV) Cluster description file. 
#               Expected columns: [Tissue, Cluster_ID, Cell_Type, Group]
# File System:  Expects a directory structure where each tissue has a subfolder 
#               containing a 'cell_data.tsv' file.
# -------------------------- #

import sys

# Hardcoded list of tissues to process, derived from a multi-line string
str_tissues = '''salivary_gland
skin
adipose_tissue
vascular
tongue
thymus
testis
stomach
spleen
small_intestine
skeletal_muscle
rectum
prostate
placenta
pbmc
pancreas
ovary
lymph_node
lung
liver
kidney
heart_muscle
fallopian_tube
eye
esophagus
endometrium
colon
bronchus
breast
brain
bone_marrow'''

l_tissues = str_tissues.split("\n")

# Path to the TSV containing biological metadata for each cluster
rna_single_cell_cluster_description = sys.argv[1]

# Dictionary to store mapped biological information
# Key: "tissue|cluster_id" -> Value: dict of metadata
d_clusters = {}

# --- STEP 1: LOAD METADATA --- #
with open(rna_single_cell_cluster_description, "r") as f:
    for line in f:
        l_line = line.strip().split("\t")
        # Normalize tissue names (lowercase and underscores) to match folder names
        tissue = l_line[0].lower().replace(" ", "_")
        cluster = l_line[1]
        
        # Build the lookup key using tissue and cluster ID
        d_clusters[tissue + "|" + cluster] = {
            "tissue": tissue, 
            "group": l_line[3], 
            "cell-type": l_line[2], 
            "cell-cluster": cluster
        }

# --- STEP 2: PROCESS TISSUE FILES --- #

# Print the header for the final integrated output
print("cell_id", "tissue", "group", "cell_type", "cluster", "umap_x", "umap_y", sep ="\t")

for tissue in l_tissues:
    # Construct path for each tissue's specific cell data (e.g., ./lung/cell_data.tsv)
    try:
        with open("./" + tissue + "/cell_data.tsv") as f:
            skip_header = f.readline() # Discard the header of the coordinate file
            
            for line in f:
                l_line = line.strip().split("\t")
                cell_id = l_line[0]
                # Format cluster ID to match the metadata keys (prefixing with 'c-')
                cluster = "c-" + l_line[1]
                umap_x = l_line[2]
                umap_y = l_line[3]
                
                # Fetch the corresponding biological annotation from our dictionary
                cell_annotation = d_clusters[tissue + "|" + cluster]
                group = cell_annotation["group"]
                cell_type = cell_annotation["cell-type"]
                
                # Output the integrated record (ID + Tissue + Biological Bio + UMAP Coordinates)
                print(cell_id, tissue, group, cell_type, cluster, umap_x, umap_y, sep = "\t")
    except FileNotFoundError:
        # Silently skip tissues if their directory or data file is missing
        continue
