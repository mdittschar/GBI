import Bio
from Bio import SeqIO
import getopt
import sys
import numpy as np
import pandas as pd

def get_seqs():
    argv = sys.argv[1:]
    print("Argv: ", argv)
    opts, args = getopt.getopt(argv, "a:b:m:s:g:", ['file1', 'file2', 'match', 'mismatch', 'gap'])
    ses = [None] * 2
    i = 0
    for opt, arg in opts:
        if opt in ["-a", "--file1", "-b", "--file2"]:
            for record in SeqIO.parse(arg, "fasta"):
                se = record.seq
                ses[i] = se
                print("Seq is: ", ses)
                i = i + 1
        elif opt == "-g" or opt =="--gap":
            gap = int(arg)
        elif opt == "-m" or opt =="--match":
            match = int(arg)
        elif opt == "-s" or opt =="--mismatch":
            mismatch = int(arg)
 
    return ses, match, mismatch, gap

def compute(sequences, match, mismatch, gap):
    rows = len(sequences[0]) + 1
    columns = len(sequences[1]) + 1
    S = np.zeros((rows, columns)).astype(int)
    T = np.zeros((rows, columns)).astype(str)
    S[0, :] = np.arange(columns)*(-gap)
    S[:, 0] = np.arange(rows)*(-gap)
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
            vleft = S[row, column-1] - gap
            vup = S[row-1, column] - gap
            values = np.array([vmatch, vleft, vup])
            vmin = np.max(values)
            S[row, column ] = vmin
            direction = directions[np.argmax(values)]
            T[row, column] = direction

    return S, T, rows, columns

    
    

def traceback(S, T, rows, columns, sequence0, sequence1):
    matrix = pd.DataFrame(S)
    tmatrix = pd.DataFrame(T)
    print(matrix)
    print(tmatrix)
    print("Sequence 0: ", sequence0)
    print("Sequence 1: ", sequence1)
    # traceback
    str_seq0 = list(str(sequence0))
    str_seq1 = list(str(sequence1))
    cur_row = rows - 1
    cur_col = columns - 1
    tstring0 = ""
    tstring1 = ""
    match_no = 0
    mismatch_no = 0
    gap_no = 0
    while cur_col != 0 and cur_row != 0:
        if T[cur_row, cur_col] == "diagonal":
            if str_seq0[cur_row -1] == str_seq1[cur_col - 1]:
                match_no = match_no + 1
            else:
                mismatch_no = mismatch_no + 1
            cur_row = cur_row - 1
            cur_col = cur_col - 1
            tstring0 = tstring0 + str_seq0[cur_row]
            tstring1 = tstring1 + str_seq1[cur_col]
        elif T[cur_row, cur_col] == "left":
            cur_col = cur_col - 1
            gap_no = gap_no + 1
            tstring0 = tstring0 + "-"
            tstring1 = tstring1 + str_seq1[cur_col]
        else:
            cur_row = cur_row - 1 
            gap_no = gap_no + 1
            tstring0 = tstring0 + str_seq0[cur_row]
            tstring1 = tstring1 + "-"
    tstring0 = tstring0[::-1]
    tstring1 = tstring1[::-1]
    print(f"Traceback strings: \n{tstring0}\n{tstring1}")

    opt_score = T[rows-1, columns-1]
    return opt_score, match_no, mismatch_no, gap_no

    

def main():
     
    sequences, match, mismatch, gap= get_seqs()
    print("Length of sequence 1: ",len(sequences[0]))
    S, T, rows, columns =  compute(sequences, match, mismatch, gap)
    opt_score, match_no, mismatch_no, gap_no = traceback(S, T, rows, columns, sequences[0], sequences[1])
    print("Optimal score: ", opt_score)
    print("Match Number: ", match_no)
    print("Mismatch Number: ", mismatch_no)
    print("Gap Number: ", gap_no)
    



if __name__ == "__main__":
    try:
        main()
    except: 
        print("Try : dittschar_auckenthaler_assignment2.py --file1 <file1> --file2 <file2> --match <match-score> --mismatch <mismatch-score> --gap <gap-score>")

    



