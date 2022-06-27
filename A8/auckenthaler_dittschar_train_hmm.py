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


allses = "".join(ses)
print(allses)