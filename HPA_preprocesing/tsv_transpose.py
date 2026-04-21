import sys

# Initialize a list to store the rows of the file
l_of_l = []

# Read the input file provided as the first command-line argument
with open(sys.argv[1], "r") as f:
    for line in f:
        # Strip newline characters and split by tabs to create a list of values
        l_line = line.strip().split("\t")
        l_of_l.append(l_line)

## CHAT GPT MATRIX TO TRANSPOSE A MATRIX STORED AS A LIST OF LISTS
def transpose_matrix(matrix):
    # Get the dimensions: num_rows is the outer list length, num_cols is the inner list length
    num_rows = len(matrix)
    num_cols = len(matrix[0])  # Assumes a rectangular matrix (consistent row lengths)
    
    # Pre-allocate the new matrix with swapped dimensions (columns become rows)
    # New shape: [num_cols][num_rows]
    transposed_matrix = [[0 for _ in range(num_rows)] for _ in range(num_cols)]
    
    # Iterate through the original rows (i) and columns (j)
    for i in range(num_rows):
        for j in range(num_cols):
            # Assign the value at [row][col] to the [col][row] position in the new matrix
            transposed_matrix[j][i] = matrix[i][j]
            
    return transposed_matrix

# Execute the transposition
transposed = transpose_matrix(l_of_l)

# Iterate through the newly created rows and print them as tab-separated strings
for l_line in transposed:
    print("\t".join(l_line))
