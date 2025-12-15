import sys
import json
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy.stats import kstest
import datetime
import os
import statistics
import math
import gc
import scanpy
import time
from datetime import datetime
import numpy as np
from scipy.sparse import csr_matrix







#-- FIXED PARAMS --#
ngenes_cutoff = 10 #Only HPOs with at least this number of related genes are going to be studied
timing = True
#------------------#


def tsv_to_csr_dense(tsv_filepath):
    """Reads a dense matrix from a TSV, skipping the first row and column (labels),
       and converts it to CSR format.

    Args:
        tsv_filepath: Path to the TSV file.

    Returns:
        A scipy.sparse.csr_matrix, or None if an error occurs.
    """
    try:
        # Use skiprows and usecols to skip the labels
        dense_matrix = np.loadtxt(tsv_filepath, delimiter='\t', skiprows=1, dtype=float)[:,1:]
        csr_matrix_result = csr_matrix(dense_matrix)
        return csr_matrix_result
    except Exception as e:
        print(f"Error reading or converting TSV: {e}")
        return None

def sparse_to_dense(sparse_list, default_value=0):
    """Converts a sparse list to a dense list.

    Args:
        sparse_list: A list of tuples, where each tuple represents a
                     non-zero element.  The first element of the tuple is
                     the index, and the second is the value.  The list
                     should be sorted by index.
        default_value: The value to use for missing elements.

    Returns:
        A dense list.
    """
    max_index = sparse_list[-1][0]  # Find the largest index
    dense_list = [default_value] * (max_index + 1)

    for index, value in sparse_list:
        dense_list[index] = value

    return dense_list






# print datetime
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print("Launched at: =", dt_string, file=sys.stderr)



HPA_file = sys.argv[1]
raw_counts = tsv_to_csr_dense(HPA_file)

with open(HPA_file, "r") as f:
	for line in f:
		l_line = line.strip().split("\t")
		l_genes_wvariants = l_line[1:]
		break

l_genes = []
for gene in l_genes_wvariants:#get rid of the variant id of the ensemble id (the part after the point)
	l_genes.append(gene.split(".")[0])


tissue_name = sys.argv[4].lower()





with open(sys.argv[3], 'r') as f:
	hpo2genes = json.load(f)


l_hpos = []
with open(sys.argv[2], "r") as f:
	for line in f:
		l_line = line.strip().split("\t")
		hpo_id = l_line[0]
		hpo_name = l_line[1]
		
		try:
			n_related_genes = len(hpo2genes[hpo_id]['EnsemblID'])
			if n_related_genes >= ngenes_cutoff:
				l_hpos.append([hpo_id, hpo_name])
		except:
			pass





print("INPUT FILES READ SUSCESFULLY!", file = sys.stderr)

print("PREPROCESSIONG...", file = sys.stderr)






## FUNCTION TO CALCULATE THE EFFECT SIZE (APLLIED TO THE KS STATISTIC, IT CORRECTS BY SAMPLE SIZE (FUENTES: LA CIBELES))

def cal_effect_size(Z, n1, n2):
	effect_size = Z*math.sqrt((n1*n2)/(n1+n2))
	return(effect_size)

n_hpos = len(l_hpos)
print("DONE!","Number of HPOs that are going to be processed: ", str(n_hpos), file=sys.stderr)


l_total_expresed = []
first_lap = True# used to count the l_total_expresed just once (as it is independent of the HPO)
n = 0
for hpoid_and_name in l_hpos:
	start_time = time.perf_counter()
	hpo = hpoid_and_name[0]#HPO identifier
	hponame = hpoid_and_name[1]#common name
	l_related_expresed = []
	n += 1
	print("Processing", hpo, "(" + hponame + ")", " [" + str(n) +"/" + str(n_hpos) + "]", end ="\n", file=sys.stderr)
	n_cell = 0
	related_genes = hpo2genes[hpo]['EnsemblID']

	first_line = True

	cell_indices = range(0, raw_counts.shape[0])
	gene_indices = range(0, raw_counts.shape[1])
	for n_cell in cell_indices:
		#print(n_cell)
		#print("Processing", hpo, "(" + hponame + "):",  "cell nº " + str(n_cell), end ="\r", file=sys.stderr) 
		n_related = 0
		n_total = 0
		non0_indices = raw_counts[n_cell].indices
		for expressed_gene_index in non0_indices:
			expressed_gene = l_genes[expressed_gene_index]
			if expressed_gene in related_genes:
				n_related += 1
		if first_lap:
			l_total_expresed.append(len(non0_indices))
		l_related_expresed.append(n_related)
	#-- CALCULATE NUMBER OF RELATED GENES TO THIS HPO --#
	n_related_genes = len(related_genes)
	
	if first_lap:
		print(f"{tissue_name}|all_genes" + "\t" + "\t".join(str(x) for x in l_total_expresed))
	print(hpo + "\t" + "\t".join(str(x) for x in l_related_expresed))
	first_lap = False
	end_time = time.perf_counter()
	elapsed_time = end_time - start_time
	print(f"Elapsed time: {elapsed_time} seconds", file=sys.stderr)
	
print("ENDED SUSCESFULLY", file=sys.stderr)
