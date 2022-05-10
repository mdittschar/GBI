import Bio
from Bio import SeqIO
import getopt
import sys
import numpy as np
import pandas as pd

def get_args():
    """
    get sequences as well as algorithm parameters from the command line

    Returns:
    --------
        ses: list of sequences
        sequences for comparison

        match: int
        match score

        mismatch: int
        mismatch score

        gap: int
        gap score
    """
    # get the command line arguments
    argv = sys.argv[1:]
    # get a tuple of the arguments to use
    opts, _ = getopt.getopt(argv, "a:b:m:s:g:", ['file1', 'file2', 'match', 'mismatch', 'gap'])
    ses = [None] * 2
    i = 0
    for opt, arg in opts:
        if opt in ["-a", "--file1", "-b", "--file2"]:
            print(f"File No {i+1}: {arg}")
            # use Biopython parser for sequences
            for record in SeqIO.parse(arg, "fasta"):
                se = record.seq
                # append sequences to list
                ses[i] = se
                i = i + 1
        # get gap, match and mismatch scores from arguments 
        elif opt == "-g" or opt =="--gap":
            print(f"Gap penalty: {arg}")
            gap = int(arg)
        elif opt == "-m" or opt =="--match":
            print(f"Match Score: {arg}")
            match = int(arg)
        elif opt == "-s" or opt =="--mismatch":
            print(f"Mismatch Score: {arg}")
            mismatch = int(arg)
 
    return ses, match, mismatch, gap

def compute(sequences, match, mismatch, gap):
    """
    Compute the sequence matrix and the traceback matrix from two sequences and linear scores

    Parameters:
    -----------
        sequences (list of sequences): sequences for alignment
        match (int): match score
        mismatch (int): mismatch score
        gap (int): gap score

    Returns:
    --------
        S (np.array): Scoring matrix for the two input sequences
        T (np.array): Traceback matrix for the two input sequences
        rows (int): number of rows (length of sequence 1 plus 1 for zero-row)
        columns (int) : number of columns (length of sequence 2 plus 1 for zero-column)
    """
    # get the number off rows and columns
    rows = len(sequences[0]) + 1
    columns = len(sequences[1]) + 1
    # initialise matrices
    S = np.zeros((rows, columns)).astype(int)
    T = np.zeros((rows, columns)).astype(str)
    S[0, :] = np.arange(columns)*(-gap)
    S[:, 0] = np.arange(rows)*(-gap)
    # direction pointers are strings
    T[0,: ] = "left"
    T[:, 0] = "up"

    # generate a list of single characters for the two sequences
    str_seq0 = list(str(sequences[0]))
    str_seq1 = list(str(sequences[1]))
    # define the possible directions for pointers
    directions = np.array(["diagonal", "left", "up"]).astype(str) # diagonal, left, up
    # all cells in the matrices are looped over
    # start from 1 because the 0th row and column are already initialised
    for column in np.arange(1,columns): 
        for row in np.arange(1,rows): 
            # get either match or mismatch score
            # since there is an additional zero-row and zero-column, reference is minus one
            if str_seq0[row- 1] == str_seq1[column- 1]:
                vmatch = S[row-1, column-1] + match
            else: 
                vmatch = S[row-1, column-1] + mismatch
            # get gap scores
            vleft = S[row, column-1] - gap
            vup = S[row-1, column] - gap
            values = np.array([vmatch, vleft, vup])
            # get the smallest score
            vmin = np.max(values)
            S[row, column] = vmin
            # point in the direciton of the lowest score
            direction = directions[np.argmax(values)]
            T[row, column] = direction

    return S, T, rows, columns

    
    

def traceback(S, T, rows, columns, sequence0, sequence1):
    """
    Trace back the optimal alignment for a given sequence and traceback matrix

    Parameters:
    -----------
        S (np.array): sequence matrix with alignment scores
        T (np.array): traceback matrix with directions
        rows (int): number of rows in S and T
        columns (int): number of columns in S and T
        sequence0 (Seq): Sequence 1 to align
        sequence1 (Seq): Sequence 2 to align

    Returns:
    --------
        opt_score (int): optimal global alignment score
        match_no (int): number of matches
        mismatch_no(int): number of mismatches 
        gap_no (int): number of gaps 
        tstring0 (str): string 1 of optimal alignment
        tstring1 (str): string 2 of optimal alignments
    """
    # matrix = pd.DataFrame(S)
    # tmatrix = pd.DataFrame(T)
    # print(matrix)
    # print(tmatrix)
    # print("Sequence 0: ", sequence0)
    # print("Sequence 1: ", sequence1)
    # traceback
    # get a characterwise list of the sequences
    str_seq0 = list(str(sequence0))
    str_seq1 = list(str(sequence1))
    # initialise current row and column
    cur_row = rows - 1
    cur_col = columns - 1
    # initialise alignment strings
    tstring0 = ""
    tstring1 = ""
    # initialise match, mismatch and gap number
    match_no = 0
    mismatch_no = 0
    gap_no = 0
    # while we have not gone through the whole matrix
    while cur_col != 0 and cur_row != 0:
        if T[cur_row, cur_col] == "diagonal":
            # if characters match, increase match number
            if str_seq0[cur_row -1] == str_seq1[cur_col - 1]:
                match_no = match_no + 1
            # else, increase mismatch number
            else:
                mismatch_no = mismatch_no + 1
            # go back up diagonally
            cur_row = cur_row - 1
            cur_col = cur_col - 1
            tstring0 = tstring0 + str_seq0[cur_row]
            tstring1 = tstring1 + str_seq1[cur_col]
        elif T[cur_row, cur_col] == "left":
            # go one cell to the left
            cur_col = cur_col - 1
            gap_no = gap_no + 1
            # add gap character to first sequence
            tstring0 = tstring0 + "-"
            tstring1 = tstring1 + str_seq1[cur_col]
        else:
            # go one cell up
            cur_row = cur_row - 1 
            gap_no = gap_no + 1
            # add gap character to second sequence
            tstring0 = tstring0 + str_seq0[cur_row]
            tstring1 = tstring1 + "-"
    # reverse the completed alignment strings
    tstring0 = tstring0[::-1]
    tstring1 = tstring1[::-1]
    
    display_seq1 = tstring1
    alignment = np.where(np.array(list(tstring0))== np.array(list(tstring1)), "|", " ")
    print(f"Traceback strings: ")
    i = 0
    while i*60 < len(tstring0):
        display_seq0 = tstring0[i*60:(i+1)*60]
        display_seq1 = tstring1[i*60:(i+1)*60]
        display_alignment = "".join(np.where(np.array(list(tstring0))== np.array(list(tstring1)), "|", " ")[i*60:(i+1)*60])
        print(f"\n{display_seq0}\n{display_alignment}\n{display_seq1}")
        i = i + 1
    
    # get the optimal alignment score
    opt_score = S[rows-1, columns-1]
    return opt_score, match_no, mismatch_no, gap_no, tstring0, tstring1

    
def visual_alignment(tstring0, tstring1):
    """
    Generate a visual alignment of two sequences following BLAST-alignment convention

    Parameters:
    -----------
        tstring0 (str): aligned string of first sequence
        tstring1 (str): aligned string of second sequence
    """
    # initialise alignment array
    alignment = np.zeros((len(tstring0), 3)). astype(str)
    # first row ist first sequence
    alignment[:,0] = list(tstring0)
    # match character where matches are present
    alignment[:,1] = np.where(np.array(list(tstring0))== np.array(list(tstring1)), "|", "")
    # second row is second sequence
    alignment[:,2] = list(tstring1)
    # generate Dataframe
    a_df = pd.DataFrame(alignment.T)
    # set visualisation parameters to fit a LARGE screen
    # according to fasta convention display 60 pairs
    # please adjust display.width if dataframe is not displayed correctly
    pd.set_option("display.max_columns", len(tstring0))
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 0)
    print(a_df)
    


def main():
     
    sequences, match, mismatch, gap= get_args()
    print("Length of sequence 1: ",len(sequences[0]))
    S, T, rows, columns =  compute(sequences, match, mismatch, gap)
    opt_score, match_no, mismatch_no, gap_no, astring0, astring1 = traceback(S, T, rows, columns, sequences[0], sequences[1])
    print("Optimal score: ", opt_score)
    print("Match Number: ", match_no)
    print("Mismatch Number: ", mismatch_no)
    print("Gap Number: ", gap_no)
    visual_alignment(astring0, astring1)

    



if __name__ == "__main__":
    try:
        main()
    except: 
        print("Try : python dittschar_auckenthaler_assignment2.py --file1 <file1> --file2 <file2> --match <match-score> --mismatch <mismatch-score> --gap <gap-score>")

    



