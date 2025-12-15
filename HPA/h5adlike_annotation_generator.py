import sys

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

rna_single_cell_cluster_description = sys.argv[1]

d_clusters = {}

with open(rna_single_cell_cluster_description, "r") as f:
	for line in f:
		l_line = line.strip().split("\t")
		tissue = l_line[0].lower().replace(" ", "_")
		cluster = l_line[1]
		d_clusters[tissue + "|" + cluster] = {"tissue": tissue, "group": l_line[3], "cell-type": l_line[2], "cell-cluster": cluster}

print("cell_id", "tissue", "group", "cell_type", "cluster", "umap_x", "umap_y", sep ="\t")
for tissue in l_tissues:
	with open("./" + tissue + "/cell_data.tsv") as f:
		skip_header = f.readline()
		for line in f:
			l_line = line.strip().split("\t")
			cell_id = l_line[0]
			cluster = "c-" + l_line[1]
			umap_x = l_line[2]
			umap_y = l_line[3]
			cell_annotation = d_clusters[tissue + "|" + cluster]
			group = cell_annotation["group"]
			cell_type = cell_annotation["cell-type"]
			print(cell_id, tissue, group, cell_type, cluster, umap_x, umap_y, sep = "\t")





