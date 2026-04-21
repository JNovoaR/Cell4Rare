# Script to extract child terms (sub-terms) of a given ontology term (.obo).
# Output is tab-separated: [ID, Name, Relationship_Type]

import argparse
import networkx
import obonet

# --- ARGUMENT PARSING --- #

parser = argparse.ArgumentParser()
parser.add_argument('obo', type=str, help='Obo file path or URL')
parser.add_argument('tissue', type=str, help='Tissue name or ID to search')
args = parser.parse_args()

# --- ONTOLOGY LOADING --- #

# Read the .obo file and convert it into a directed graph
graph = obonet.read_obo(args.obo) 

# Create lookup dictionaries for ID-to-Name and Name-to-ID conversion
id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True) if 'name' in data}
node_data = graph.nodes(data=True)

# --- TERM IDENTIFICATION --- #

# Attempt to find the specific ID if a common name was provided as input
try:
    tissue = name_to_id[args.tissue]
except KeyError:
    tissue = args.tissue # Use as-is if already an ID

# --- RECURSIVE SEARCH FUNCTION --- #

def get_file(tissue, keylist=None, exclude="ignore"):
    """
    Traverses the ontology graph to find all descendants of a tissue.
    
    Args:
        tissue: Starting node ID.
        keylist: List of edge types to follow (e.g., 'is_a', 'part_of').
        exclude: Behavior for filtering edge types.
    """
    id_list = [tissue]
    print(f'{tissue}\t{id_to_name[tissue]}')
    
    used_ids = []      # Track nodes already processed
    temp_list = []     # Store children found in the current iteration
    repeated_list = [] # Prevent infinite loops/redundancy
    terms_avaiable = True

    # Breath-first search traversal
    while terms_avaiable:
        for ub_id in id_list:
            if ub_id not in used_ids:
                # Explore edges coming into the node (these are the 'children' in obonet)
                for parent, child, key in graph.in_edges(ub_id, keys=True):
                    # Filter based on relationship keys (e.g., only follow 'is_a')
                    if exclude == "ignore": 
                        condition = True
                    elif exclude == True: 
                        condition = key not in keylist
                    elif exclude == False: 
                        condition = key in keylist
                    
                    if condition:
                        if parent not in repeated_list:
                            # Print the found child term and its relationship
                            print(f'{parent}\t{id_to_name[parent]}\t{key}')
                            temp_list.append(parent)
                            repeated_list.append(parent)
                
                used_ids.append(ub_id)
        
        # Move to the next "generation" of children
        id_list = temp_list
        temp_list = []
        
        # Stop if no more children are found
        if len(id_list) == 0: 
            terms_avaiable = False

# --- EXECUTION --- #

# Define which biological relationships are valid to follow
# 'is_a' follows hierarchy (e.g., Lung is_a Organ)
# 'part_of' follows anatomy (e.g., Alveolus part_of Lung)
keys = ('is_a', 'part_of')

# Run the extraction starting from the target tissue
get_file(tissue, keys, exclude=False)
