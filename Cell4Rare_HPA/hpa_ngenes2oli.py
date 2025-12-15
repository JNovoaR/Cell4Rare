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







def convert_to_integers(string_list):
    try:
        integer_list = [int(item) for item in string_list]
        return integer_list
    except ValueError as e:
        print(f"Error converting list to integers: {e}")
        return None

## FUNCTION TO CALCULATE THE EFFECT SIZE (APLLIED TO THE KS STATISTIC, IT CORRECTS BY SAMPLE SIZE (FUENTES: LA CIBELES))

def cal_effect_size(Z, n1, n2):
	effect_size = Z/math.sqrt((n1*n2)/(n1+n2))
	return(effect_size)

#-- FIXED PARAMS --#
ngenes_cutoff = 10 #Only HPOs with at least this number of related genes are going to be studied

#------------------#



f_cell_annotation = sys.argv[1] #h5adlike_annotation_generator.py >hpa_annotation.tsv

f_rel_hpos = sys.argv[2]




l_hpos = []
with open(f_rel_hpos, "r") as f:
	for line in f:
		l_line = line.strip().split("\t")
		hpo_id = l_line[0]
		hpo_name = l_line[1]
		
		
		n_related_genes = int(l_line[3])
		if n_related_genes >= ngenes_cutoff:
			l_hpos.append([hpo_id, hpo_name, n_related_genes])

l_tissues = []
l_total_expresed = []#background
for fname in sys.argv[3:]:
	with open(fname, "r") as f:
		line = f.readline()
		l_line = line.strip().split("\t")
		tissue = l_line[0].split("|")[0].lower().replace(" ", "_")
		l_tissues.append(tissue)
		l_total_expresed += convert_to_integers(l_line[1:])


n_cell = 0
d_tissue2cell = {}
d_groupntissue2cell = {}
d_celltypengroupntissue2cell = {}
d_clusterncelltypengroupntissue2cell = {}


for included_tissue in l_tissues:
	with open(f_cell_annotation, "r") as f:
		for line in f:
			if "\"" in line:
				new_line = ""
				inside_quotation = False
				for char in line:
					if char == "\"":
						char = ""
						if inside_quotation:
							inside_quotation = False
						else:
							inside_quotation = True
					if char == ",":
						if inside_quotation:
							char = ";"
					new_line += char
				line = new_line
			### RELEVANT TO ARODRIGUEZ ###
			l_line = line.strip().split("\t")
			cell_id = l_line[0].lower()
			tissue = l_line[1].lower()
			group = l_line[2].lower()
			cell_type = l_line[3].lower()
			cluster = l_line[4].lower()
			if tissue == included_tissue:
				if tissue not in d_tissue2cell.keys():
					d_tissue2cell[tissue] = [n_cell]
				else:
					d_tissue2cell[tissue].append(n_cell)
					
				if f"{tissue}/{group}" not in d_groupntissue2cell.keys():
					d_groupntissue2cell[f"{tissue}/{group}"] = [n_cell]
				else:
					d_groupntissue2cell[f"{tissue}/{group}"].append(n_cell)
					
				if f"{tissue}/{group}/{cell_type}" not in d_celltypengroupntissue2cell.keys():
					d_celltypengroupntissue2cell[f"{tissue}/{group}/{cell_type}"] = [n_cell]
				else:
					d_celltypengroupntissue2cell[f"{tissue}/{group}/{cell_type}"].append(n_cell)
				
				if f"{tissue}/{group}/{cell_type}/{cluster}" not in d_clusterncelltypengroupntissue2cell.keys():
					d_clusterncelltypengroupntissue2cell[f"{tissue}/{group}/{cell_type}/{cluster}"] = [n_cell]
				else:
					d_clusterncelltypengroupntissue2cell[f"{tissue}/{group}/{cell_type}/{cluster}"].append(n_cell)
				
				n_cell += 1

			### ---------------------- ###

n_hpos = len(l_hpos)
n = 1

#print(d_groupntissue2cell["salivary_gland/glandular epithelial cells"])

for id_n_name in l_hpos:
	hpo_id = id_n_name[0]
	hpo_name = id_n_name[1]
	n_related_genes = id_n_name[2]
	l_related_expresed = []#distribution
	for fname, tissue in zip(sys.argv[3:], l_tissues):#loop through ngenes files and its corresponding tissues in paralel
		with open(fname, "r") as f:
			hpo_ispresent = False
			for line in f:
				if line.startswith(hpo_id):
					l_line = line.strip().split("\t")
					l_related_expresed += convert_to_integers(l_line[1:])
					hpo_ispresent = True
					break
			if not hpo_ispresent:
				print(f"ERROR: at least one of the HPOs ({hpo_id}:{hpo_name} present in the input ('relHPO_tissue') is not present in at least one of the ngenes files. Exiting...", file = sys.stderr)
				exit()
	
	print("Processing", hpo_id, "(" + hpo_name + ")", " [" + str(n) +"/" + str(n_hpos) + "]", end ="\n", file=sys.stderr)
	##-- Check vectors sizes --#
	if len(l_related_expresed) != len(l_total_expresed):
		print("ERROR: l_related_expresed is not equal in lenght to l_total_expresed. I cant think of a situation where this would happend. :(", file = sys.stderr)
		exit()
	
	
	#-- FIT A LINEAR REGRESSION --#
	X_tofit = np.array(l_total_expresed).reshape(-1, 1)
	y_tofit = np.array(l_related_expresed)
	model = LinearRegression()
	model.fit(X_tofit, y_tofit)
	m = model.coef_[0]
	b = model.intercept_
	#-- SCATTER PLOT WITH LINEAR REGRESSION --#
	#x_line = np.linspace(min(X), max(X), 100)  # Generate x values for the line
	#y_line = m * x_line + b  # Calculate y values using the equation
	#save_scatter_plot(l_total_expresed, l_related_expresed, "Nº total genes expressed", hpo, x_line, y_line,  "Nº related genes expressed", "prueba_" + hpo + ".jpg")
	
	#-- CALCULATE DISTANCE BETWEEN OBSERVED Y AND EXPECTED Y VALUE IN THE LINEAR REGRESION (CORRESPONDING TO THE SAME X) --#
	distances = []
	for x, y_obs in zip(l_total_expresed, l_related_expresed):
		y_exp = x*m + b
		dist = y_obs - y_exp
		distances.append(dist)
	### RELEVANT TO ARODRIGUEZ ###
	list_of_dictionaries = [d_tissue2cell, d_groupntissue2cell, d_celltypengroupntissue2cell, d_clusterncelltypengroupntissue2cell]
	list_of_agrupation_levels = ["tissue", "group", "cell-type", "cluster"]
	### ---------------------- ###
	for d_agrupation2cell, agrupation_level in zip(list_of_dictionaries, list_of_agrupation_levels):
		for agrupation in d_agrupation2cell.keys():
			agrupation_indices = d_agrupation2cell[agrupation]
			agrupation_distances = [distances[i] for i in range(len(distances)) if i in agrupation_indices]
			rest_distances = [distances[i] for i in range(len(distances)) if i not in agrupation_indices]
			if len(rest_distances) == 0: #in casses in wich there is only one tissue studied, or only one group in the whole tissue and only a tissue studie (the agrupation == all cells)
				ks_statistic = "NA"
				p_value = "NA"
				significative = "NA"
				effect_size = "NA"
				
			else:
				ks_statistic, p_value = kstest(agrupation_distances, rest_distances, alternative = 'less')#We ask: is sample greater than background?
				significative = int(p_value < 0.01)
				effect_size = cal_effect_size(ks_statistic, len(agrupation_distances), len(rest_distances))
			#print(sorted(distances) == sorted(rest_distances + cluster_distances))
			#print(sorted(cluster_indices))
			print(hpo_id, hpo_name, n_related_genes, agrupation_level, agrupation, len(agrupation_distances), len(rest_distances), ks_statistic, p_value, effect_size, significative, sep="\t")
	n += 1
		
				
	
	
	
