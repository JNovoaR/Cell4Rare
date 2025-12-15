import sys

l_of_l = []
with open(sys.argv[1], "r") as f:
	for line in f:
		l_line = line.strip().split("\t")
		l_of_l.append(l_line)

##CHAT GPT MATRIX TO TRASPOSE A MATRIX STORED AS A LIST OF LISTS
def transpose_matrix(matrix):
	# Get the number of rows and columns in the original matrix
	num_rows = len(matrix)
	num_cols = len(matrix[0])  # Assuming all rows have the same number of columns
	# Create a new matrix to store the transposed matrix
	transposed_matrix = [[0 for _ in range(num_rows)] for _ in range(num_cols)]
	# Iterate through the original matrix and fill the transposed matrix
	for i in range(num_rows):
		for j in range(num_cols):
			transposed_matrix[j][i] = matrix[i][j]
	return transposed_matrix

transposed = transpose_matrix(l_of_l)

for l_line in transposed:
	print("\t".join(l_line))
