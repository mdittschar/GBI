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
    S = np.zeros((rows, columns)).astype(int)
    T = np.zeros((rows, columns)).astype(str)
    S[0, :] = np.arange(columns)*gap
    S[:, 0] = np.arange(rows)*gap
    T[0,: ] = "left"
    T[:, 0] = "up"

    
    str_seq0 = list(str(sequences[0]))
    str_seq1 = list(str(sequences[1]))
    no_match = 0
    no_mismatch = 0
    directions = np.array(["diagonal", "left", "up"]).astype(str) # diagonal, left, up
    for column in np.arange(1,columns): 
        for row in np.arange(1,rows): 
            if str_seq0[row- 1] == str_seq1[column- 1]:
                no_match = no_match + 1
                vmatch = S[row-1, column-1] + match
            else: 
                vmatch = S[row-1, column-1] + mismatch
                no_mismatch = no_mismatch + 1
            vleft = S[row, column-1] + gap
            vup = S[row-1, column] + gap
            values = np.array([vmatch, vleft, vup])
            vmin = np.min(values)
            S[row, column ] = vmin
            direction = directions[np.argmin(values)]
            T[row, column] = direction

    
    matrix = pd.DataFrame(S)
    tmatrix = pd.DataFrame(T)
    print(matrix)
    print(tmatrix)
    print("Sequence 0: ", sequences[0])
    print("Sequence 1: ", sequences[1])
    # traceback
    
    cur_row = rows - 1
    cur_col = columns - 1
    tstring0 = ""
    tstring1 = ""
    while cur_col != 0 and cur_row != 0:
        if T[cur_row, cur_col] == "diagonal":
            cur_row = cur_row - 1
            cur_col = cur_col - 1
            tstring0 = tstring0 + str_seq0[cur_row]
            tstring1 = tstring1 + str_seq1[cur_col]
        elif T[cur_row, cur_col] == "left":
            cur_col = cur_col - 1
            tstring0 = tstring0 + "-"
            tstring1 = tstring1 + str_seq1[cur_col]
        else:
            cur_row = cur_row - 1 
            tstring0 = tstring0 + str_seq0[cur_row]
            tstring1 = tstring1 + "-"
    tstring0 = tstring0[::-1]
    tstring1 = tstring1[::-1]
    print(f"Traceback strings: \n{tstring0}\n{tstring1}")





if __name__ == "__main__":
    
    sequences = get_seqs()
    print("Length of sequence 1: ",len(sequences[0]))
    compute(sequences, -2, 2, 4)

# Needleman-Wunsch

