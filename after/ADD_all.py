import sys
import os
import numpy as np
from statsmodels.stats.multitest import multipletests
import networkx
import obonet

f_oli = sys.argv[1]
folder_ngenes = sys.argv[2]
ts_or_hpa = sys.argv[3]

hpo_col = "NA"
if ts_or_hpa == "hpa":
	hpo_col = 1
	f_cell2coment = "../validation/hpa2coment"
elif ts_or_hpa == "ts":
	hpo_col = 0
	f_cell2coment = "../validation/ts2coment"
else:
	print("Error: argv[3] must be 'hpa' or 'ts'", file=sys.stderr)
	exit()

maxngenes_threshold = 10
sensitive_2_nonexistant_in_coment = True #true to consider as "NA" the row in wich the HPO has 0 matches wich any CL
f_raw_coment = "../../coment_net/DATA/ALL_CL-HPO_s"
hpo_graph = obonet.read_obo("../DATA/HPO/hp.obo") #read hpobo
topparent_term = "HP:0000118"#(Phenotypic abnormality)


ngenes_col = hpo_col + 2
cluster_lvl_col = hpo_col + 3
cluster_col = hpo_col + 4
pval_col = hpo_col + 8

###l_tissues (to know the names of the ngenes of the files)
l_tissues = []
with open(f_oli, "r") as f:
	for line in f:
		l_line = line.strip().split("\t")
		if ts_or_hpa == "hpa":
			tissue = l_line[hpo_col - 1]
		elif ts_or_hpa == "ts":
			tissue = l_line[cluster_col].split("/")[0]
		if tissue not in l_tissues:
			l_tissues.append(tissue)


###--- MAXGENES ---###
d_hpo2ngenes = {}
for tissue in l_tissues:
	d_hpo2ngenes[tissue] = {}
	if ts_or_hpa == "ts":
		ngenes_filename = "/ngenes_" + tissue.title().replace(" ", "_")
	else:
		ngenes_filename = "/ngenes_" + tissue.lower().replace(" ", "_")
	with open(folder_ngenes + ngenes_filename, "r") as f:
		for line in f:
			l_line = line.strip().split("\t")
			hpo = l_line[0]
			max_genes = max([int(x) for x in l_line[1:]])
			d_hpo2ngenes[tissue][hpo] = max_genes

l_f_mg = []
with open(f_oli, "r") as f:
	for line in f:
		l_line = line.strip().split("\t")
		hpo = l_line[hpo_col]
		if ts_or_hpa == "hpa":
			tissue = l_line[hpo_col - 1]
		elif ts_or_hpa == "ts":
			tissue = l_line[cluster_col].split("/")[0]
		if d_hpo2ngenes[tissue][hpo] >= maxngenes_threshold:
			max_ngenes = d_hpo2ngenes[tissue][hpo]
			l_f_mg.append("\t".join(l_line[:ngenes_col]) +"\t"+ str(d_hpo2ngenes[tissue][hpo]) +"\t"+ "\t".join(l_line[ngenes_col + 1:]))


#for line in l_f_mg:
#	print(line)

###--- FDR CORRECTION ---###

d_pvals = {}
d_counter = {}
if ts_or_hpa == "hpa":
	cluster_levels = ["type-group", "cell-type", "cell-cluster"]
	
	for i in cluster_levels:
		d_pvals[i] = []
		d_counter[i] = 0
elif ts_or_hpa == "ts":
	cluster_levels = ["group", "cell-type"]
	for i in cluster_levels:
		d_pvals[i] = []
		d_counter[i] = 0
else:
	print("Error: arg 3 must be 'hpa' or 'ts'", file = sys.stderr)
	exit()



for line in l_f_mg:
	l_line = line.strip().split("\t")
	pval = l_line[pval_col]
	cluster_lvl = l_line[cluster_lvl_col]
	if pval != "NA":
		pval = float(pval)
		d_pvals[cluster_lvl].append(pval)

# Apply FDR correction
d_pvals_corrected = {}

for cluster_lvl in cluster_levels:
	l_reject, d_pvals_corrected[cluster_lvl], _, _ = multipletests(d_pvals[cluster_lvl], alpha=0.05, method='fdr_bh')


l_f_mg_fdr = []


for line in l_f_mg:
	l_line = line.strip().split("\t")
	pval = l_line[pval_col]
	cluster_lvl = l_line[cluster_lvl_col]
	if pval != "NA":
		l_f_mg_fdr.append(line.strip() +"\t"+ str(d_pvals_corrected[cluster_lvl][d_counter[cluster_lvl]]))
		d_counter[cluster_lvl]+=1
	else:
		l_f_mg_fdr.append(line.strip() +"\t"+ "NA")

#for line in l_f_mg_fdr:
#	print(line)

###--- COMENT TAG ---###



d_hpa2coment = {}

l_cl_present = []
l_hpo_present = []

with open(f_raw_coment, "r") as f:
	for line in f:
		l_line = line.strip().split("\t")
		l_cl_present.append(l_line[0])
		l_hpo_present.append(l_line[1])

with open(f_cell2coment, "r") as f:
	for line in f:
		l_line = line.strip().split("\t")
		if l_line[0].lower().strip() not in d_hpa2coment.keys():
			d_hpa2coment[l_line[0].lower().strip()] = [l_line[2]]
		else:
			d_hpa2coment[l_line[0].lower().strip()].append(l_line[2])



l_f_mg_fdr_cmnt = []


for line in l_f_mg_fdr:
	coment_result = "Error!"
	l_line = line.strip().split("\t")
	if l_line[hpo_col] not in l_hpo_present and sensitive_2_nonexistant_in_coment:
		coment_result = "NA"
	else:
		match = False
		exist_in_coment = True
		if l_line[cluster_col - 1] == "cell-cluster":
			cell_type = l_line[cluster_col].split("(")[1].split(")")[0]
			if cell_type in d_hpa2coment.keys():
				if l_line[hpo_col] in d_hpa2coment[cell_type]:
					coment_result = "CoMent"
				else:
					coment_result = "-"
			else:
				coment_result = "NA"
				
		elif l_line[cluster_col - 1] == "cell-type":
			if "/" in l_line[cluster_col]:
				cell_type = l_line[cluster_col].split("/")[-1]
			else:
				cell_type = l_line[cluster_col]
			if cell_type in d_hpa2coment.keys():
				if l_line[hpo_col] in d_hpa2coment[cell_type]:
					coment_result = "CoMent"
				else:
					coment_result = "-"
			else:
				coment_result = "NA"
		elif l_line[cluster_col - 1] in ["type-group", "tissue","group"]:
			coment_result = "NA"
		else:
			print("error", file = sys.stderr)
			exit()
	l_f_mg_fdr_cmnt.append(line.strip() + "\t" + coment_result)
			
print(f"Fixed param: sensitive_2_nonexistant_in_coment = {sensitive_2_nonexistant_in_coment}", file = sys.stderr)


###--- SPECIFICITY ---###

l_f_mg_fdr_cmnt_spec = []


def get_children(graph, term):
	id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data = True)}

	l_terms = [term]


	for child, parent, key in graph.in_edges(term, keys=True):
		if child not in l_terms:
			l_terms.append(child)

	return(l_terms)

### The code was changed to find the specificity of an hpo within a tissue, since due to the fdr I need to run
### this code over the whole concatenation of all the tissue files, not tissue by tissue.
 


d_hpos_intissue = {}	

for line in l_f_mg_fdr_cmnt:
	l_line = line.strip().split("\t")
	if ts_or_hpa == "hpa":
		tissue = l_line[hpo_col - 1]
	elif ts_or_hpa == "ts":
		tissue = l_line[cluster_col].split("/")[0]
	if tissue in d_hpos_intissue:
		d_hpos_intissue[tissue].append(l_line[hpo_col])
	else:
		d_hpos_intissue[tissue] = [l_line[hpo_col]]


for key in d_hpos_intissue:
	d_hpos_intissue[key] = list(set(d_hpos_intissue[key]))


for line in l_f_mg_fdr_cmnt:
	l_line = line.strip().split("\t")
	if ts_or_hpa == "hpa":
		tissue = l_line[hpo_col - 1]
	elif ts_or_hpa == "ts":
		tissue = l_line[cluster_col].split("/")[0]
	hpo_id = l_line[hpo_col]
	l_children = get_children(hpo_graph, hpo_id)
	intersec = list(set(l_children) & set(d_hpos_intissue[tissue]))
	spec = len(intersec)-1#the parent is included in the list of children, so we sustract 1
	l_f_mg_fdr_cmnt_spec.append(line.strip() + "\t" + str(spec))


###--- TOPTERM ---###

id_to_name = {id_: data.get('name') for id_, data in hpo_graph.nodes(data = True)}

def get_children2(graph, term):
	id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data = True)}

	l_terms = [term]

	for term in l_terms:
		for child, parent, key in graph.in_edges(term, keys=True):
			if child not in l_terms:
				l_terms.append(child)

	return(l_terms)


top_children = []
for child, parent, key in hpo_graph.in_edges(topparent_term, keys=True):
	top_children.append([child, id_to_name[child]])


prev_hpo = ""
for line in l_f_mg_fdr_cmnt_spec:
	l_line = line.strip().split("\t")
	hpo = l_line[hpo_col]
	if hpo != prev_hpo:
		prev_hpo = hpo
		corresponding_topterms = []
		for topterm in top_children:
			children = get_children2(hpo_graph, topterm[0])
			if hpo in children:
				corresponding_topterms.append(topterm[1])

	print(line.strip() + "\t" + "\t".join(corresponding_topterms))


