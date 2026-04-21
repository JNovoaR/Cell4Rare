# --- INPUT DESCRIPTION --- #
# sys.argv[1]: (TSV) List of anatomy/cell terms (e.g., UBERON IDs).
# sys.argv[2]: (OBO) The HPO ontology file (hp.obo).
# sys.argv[3]: (JSON) Mapping file from Anatomy terms to HPO terms (uberon2hpo).
# sys.argv[4]: (JSON) Mapping file from HPO terms to Gene IDs (hpo2genes).
# sys.argv[5]: (TSV) Manual HPO annotations for tissues.
# sys.argv[6]: (String) The name of the specific tissue to process.
# -------------------------- #

import sys
import json
import datetime
import networkx
import obonet

# --- LOAD INPUTS --- #

# 1. Load the anatomy/child terms from the first argument
l_children = []
with open(sys.argv[1], "r") as f:
    for line in f:
        l_line = line.strip().split("\t")
        l_children.append(l_line[0])

# 2. Load the HPO graph structure
graph = obonet.read_obo(sys.argv[2]) 

# 3. Load the Anatomy-to-HPO cross-reference dictionary
with open(sys.argv[3], 'r') as f:
    ub2hpo = json.load(f)

tissue_name = sys.argv[6]
print("Processing " + tissue_name, file=sys.stderr)

# 4. Load Manual HPO Annotations
# Filters the manual file to find HPO terms specifically labeled for the current tissue
l_manual_HPO = []
with open(sys.argv[5], 'r') as f:
    for line in f:
        if line.startswith("#") or line.strip() == "":
            continue
        l_line = line.strip().split("\t")
        hpo_id = l_line[1]
        
        # Check if tissue name matches (handling spaces vs underscores)
        if l_line[0] == tissue_name:
            l_manual_HPO.append(hpo_id)
        else:
            if l_line[0].split(" ") == tissue_name.split("_") or l_line[0].split(" ") == tissue_name.split("-"):
                l_manual_HPO.append(hpo_id)

# 5. Load HPO-to-Genes mapping
with open(sys.argv[4], 'r') as f:
    hpo2genes = json.load(f)

# --- PROCESS AUTOMATED MAPPINGS --- #

# Find HPO terms related to the tissue via the UBERON/anatomy IDs
l_annotated_HPO = []
for ub in l_children:
    try:
        hpos = ub2hpo[ub]
    except KeyError:
        hpos = []
    for hpo in hpos:
        l_annotated_HPO.append(hpo)

# --- ONTOLOGY EXPANSION & FILTERING --- #

def the_rest(l_input):
    """
    Expands a list of HPO terms to include all descendants (children) 
    and validates them against available gene data.
    """
    l_to_print = []
    l_nogenes = []
    
    # RECURSIVE EXPANSION:
    # Uses the graph to find all "in_edges" (children in OBO format)
    # This ensures that if 'Heart Disease' is mapped, all sub-types are also included.
    l_related_hpos_wchildren = l_input
    for term in l_related_hpos_wchildren:
        for child, parent, key in graph.in_edges(term, keys=True):
            if child not in l_related_hpos_wchildren:
                l_related_hpos_wchildren.append(child)

    # GENE VALIDATION:
    # Checks if the expanded HPO terms have associated Ensembl IDs.
    for hpo in l_related_hpos_wchildren:
        try:
            related_genes = hpo2genes[hpo]['EnsemblID']
            hpo_name = hpo2genes[hpo]['HPOname']
            l_to_print.append([hpo, hpo_name, str(len(related_genes))])    
        except KeyError:
            # If no genes are related, the HPO is discarded
            l_nogenes.append(hpo)
            
    return l_to_print, l_nogenes

# Run expansion for both automated (UBERON) and manual lists
l_annotated_HPO_wchildren, l1_nogenes = the_rest(l_annotated_HPO)
l_manual_HPO_wchildren, l2_nogenes = the_rest(l_manual_HPO)

# Count total unique HPOs discarded due to lack of genetic evidence
n_nogenes = len(set(l1_nogenes + l2_nogenes))

# --- FINAL OUTPUT --- #

printed = []

# Print terms found via automated UBERON mapping
for i in l_annotated_HPO_wchildren:
    if i in l_manual_HPO_wchildren:
        # Term exists in both automated and manual lists
        print(i[0], i[1], "BOTH", i[2], sep="\t")
        printed.append(i)
    else:
        # Term found only via UBERON
        print(i[0], i[1], "UBERON", i[2], sep="\t")

# Print remaining terms found only via Manual mapping
for i in l_manual_HPO_wchildren:
    if i not in printed:
        print(i[0], i[1], "MANUAL", i[2], sep="\t")

print("Done!", n_nogenes, "HPOs were discarded because they had no related genes\n", file=sys.stderr)
