import numpy as np
import sys
import getopt
import Bio
from Bio import SeqIO


# fasta reader from console taken from our solution for A3
argv = sys.argv[1:]
# get a tuple of the arguments to use
opts, _ = getopt.getopt(argv, "f:", ['file'])
ses = [None] * 10

ids = []
i = 0
for opt, arg in opts:
    if opt in ["-f", "--file1"]:
        print(f"File No {i+1}: {arg}")
        # use Biopython parser for sequences
        for record in SeqIO.parse(arg, "fasta"):
            print("record: ", i)
            se = record.seq
            
            id = record.id
            ids = np.append(ids, id)

            ses[i] = str(se)

            i = i + 1


allses = np.array(list("".join(ses)))

trans_mat = np.zeros((2,2))
for s in ses:
    for i, n in enumerate(s):
        if len(s) - 1 > i:
            if n == "G" and s[i+1] == "G":
                trans_mat[0,0] = trans_mat[0,0] + 1
            elif n == "C" and s[i+1] == "G":
                trans_mat[1,0] = trans_mat[1,0] + 1
            elif n == "C" and s[i+1] == "C":
                trans_mat[1,1] = trans_mat[1,1] + 1
            else:
                trans_mat[0,1] = trans_mat[0,1] + 1

trans_mat[0,:] = trans_mat[0,:]/sum(trans_mat[0,:])
trans_mat[1,:] = trans_mat[1,:]/sum(trans_mat[1,:])
print(trans_mat)