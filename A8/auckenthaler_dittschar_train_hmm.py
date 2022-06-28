import numpy as np
import sys
import getopt
import Bio
from Bio import SeqIO


# fasta reader from console taken from our solution for A3
argv = sys.argv[1:]
# get a tuple of the arguments to use
opts, _ = getopt.getopt(argv, "f:n:", ['file', 'name'])
ses = [None] * 10

i = 0
for opt, arg in opts:
    if opt in ["-f", "--file1"]:
        # use Biopython parser for sequences
        for record in SeqIO.parse(arg, "fasta"):
            se = record.seq
            # get sequences
            ses[i] = str(se)

            i = i + 1
    # get new transition matrix name
    elif opt in ["-n", "--name"]:
        matname = arg
# # join all sequences together for probability computation
# allses = np.array(list("".join(ses)))

# define number of rows
rows = 4
trans_mat = np.zeros((rows, rows))
# for every sequence
for s in ses:
    # for every character in each sequence
    for i, n in enumerate(s):
        # if it's the first character, append to "b" row
        if i == 0:
            if n == "G":
                trans_mat[2,0] = trans_mat[2,0] + 1
            else: 
                trans_mat[2,1] = trans_mat[2,1] + 1
        # if it's not the end character, append to corresponding row and column
        elif len(s) - 1 > i:
            if n == "G" and s[i+1] == "G":
                trans_mat[0,0] = trans_mat[0,0] + 1
            elif n == "C" and s[i+1] == "G":
                trans_mat[1,0] = trans_mat[1,0] + 1
            elif n == "C" and s[i+1] == "C":
                trans_mat[1,1] = trans_mat[1,1] + 1
            else:
                trans_mat[0,1] = trans_mat[0,1] + 1
        # in the end, append transition from nucleotide to end
        else:
            if n == "G":
                trans_mat[0, 3] = trans_mat[0,3] + 1
            else: 
                trans_mat[1, 3] = trans_mat[1, 3] + 1
            trans_mat[3, 3] = 1
# divide by the total of each row
for row in np.arange(rows):
    trans_mat[row,:] = trans_mat[row,:]/sum(trans_mat[row,:])
# # replace cells in transition matrix by very small value where it is zero
# trans_mat = np.where(trans_mat == 0, np.finfo(float).eps, trans_mat)
# format the output string
output_string = f"# The transition matrix can be found in the numpy file: {matname}_only.npy\n\
# Number of states:\n\
{rows}\n\
# State labels: *=b, +=e\n\
G C * + \n\
# Transition matrix P:\n\
{trans_mat}"

# write to file
with open(f"{matname}.txt", "w") as f: 
    f.write(output_string)
print(output_string)
np.save(f"{matname}_only", trans_mat)