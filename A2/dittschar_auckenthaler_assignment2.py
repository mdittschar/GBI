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

    # get the optimal alignment score
    opt_score = S[rows-1, columns-1]
    return opt_score, match_no, mismatch_no, gap_no, tstring0, tstring1

    
def visual_alignment(tstring0, tstring1,filename, match_no, mismatch_no, gap_no, match, mismatch, gap, pair_no=60 ):
    """
    Generate a visual alignment of two sequences following BLAST-alignment convention and write it to a text file
    (TASK 5)

    Parameters:
    -----------
        tstring0 (str): aligned string of first sequence
        tstring1 (str): aligned string of second sequence
        pair_no (int): number of aligned pairs to show in one row
        filename (String): tet.file name 
        match_no, (int): no. of matches in alignment
        mismatch_no (int): no of mismatches in alignment
        gap_no (int):   no. of gaps in alignment 
        match (int): match-score parameter
        mismatch (int): mismatch- score parameter
        gap (int): gap-score parameter
    """
    print(f"Traceback strings: ")
    i = 0
    seq_lens = len(tstring0)
    # continue while alignment is not fully visualised yet

    with open(filename, 'w') as file_out:

        while i*pair_no < seq_lens:
            # assign strings to be visualised
            # pair_no is number of pairs to visualise in one row
            display_seq0 = tstring0[i*pair_no:(i+1)*pair_no]
            display_seq1 = tstring1[i*pair_no:(i+1)*pair_no]
            # show alignment lines where matches are present
            display_alignment = "".join(np.where(np.array(list(tstring0))== np.array(list(tstring1)), "|", " ")[i*pair_no:(i+1)*pair_no])
            # visualisation for full lines
            if i*pair_no < seq_lens - pair_no:
                empty = (pair_no - 2 - i)* " "
                # numbers at the ends of the lines to show number of pairs
                alignment_nos = "".join([str(i*pair_no + 1), empty, str((i+1)*pair_no)])
            else:
                # visualisation if there's an incomplete line left
                empty = " "*(seq_lens - pair_no*i)
                alignment_nos = "".join([str(i*pair_no + 1), empty, str(seq_lens)])
            # print the combined information
            print(f"\n{alignment_nos}\n{display_seq0}\n{display_alignment}\n{display_seq1}")
            #Task 5 a) write it to textfile
            file_out.write(f"\n{alignment_nos}\n{display_seq0}\n{display_alignment}\n{display_seq1}")
        
            i = i + 1
        #Task 5 b) and c) write parameter to text file:
        print("Match Score: ",match)
        print("Mismatch Score: ", mismatch)
        print("Gap penalty: ",gap )
        print("Match No.: ", match_no)
        print("Mismatch No.: ", mismatch_no)
        print("Gap No.: ", gap_no) 
        file_out.write(f"\nMatch Score: {match}\nMismatch Score: {mismatch}\nGap Penalty: {gap}\nMatch No.: {match_no}\nMismatch No.:{mismatch_no}\nGap No.: {gap_no}")
        file_out.close
      

def main():
     
    sequences, match, mismatch, gap= get_args()
    # get matrices and number of rows/columns
    S, T, rows, columns =  compute(sequences, match, mismatch, gap)
    # get optimal alignment score as well as number of matches, mismatches and gaps and aligned strings
    opt_score, match_no, mismatch_no, gap_no, astring0, astring1 = traceback(S, T, rows, columns, sequences[0], sequences[1])
    print("Optimal score: ", opt_score)
    #call function to print and write output of needleman-wunsch
    visual_alignment(astring0, astring1, "dittschar_auckenthaler_assignment2_global_alignment.txt", match_no, mismatch_no, gap_no, match, mismatch, gap)

    



if __name__ == "__main__":
    try:
        main()
    except: 
        print("Try : python dittschar_auckenthaler_assignment2.py -a <file1> -b <file2> -m <match-score> -s <mismatch-score> -g <gap-score>")

    



