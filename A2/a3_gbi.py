import Bio
from Bio import SeqIO
import getopt
import sys
import numpy as np
import pandas as pd

def get_seqs():
    ses = [None] * (len(sys.argv)-1)
    for i, val in enumerate(sys.argv[1:]):
        for record in SeqIO.parse(val, "fasta"):
            se = record.seq

            ses[i] = se
   
    return ses

def compute(sequences, match, mismatch, gap):
    rows = len(sequences[0]) + 1
    columns = len(sequences[1]) + 1
    str_index = list(str(sequences[0]))
    str_index.insert(0,'0')
    str_columns = list(str(sequences[1]))
    str_columns.insert(0, '0')
    matrix = pd.DataFrame(columns = str_columns, index=str_index)
    print(matrix)

    matrix.loc["0"] = np.arange(len(str_columns)) * gap
    matrix.loc[:,"0"] = np.arange(len(str_index)) * gap
    print(matrix)


if __name__ == "__main__":
    
    sequences = get_seqs()
    print("Length of sequence 1: ",len(sequences[0]))
    compute(sequences, -2, 2, 4)

# Needleman-Wunsch

