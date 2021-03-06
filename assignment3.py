import Bio
from Bio import SeqIO
import getopt
import sys
import numpy as np
import pandas as pd
from itertools import combinations, product
from operator import itemgetter

def get_args():
    """
    get sequences as well as algorithm parameters from the command line

    Returns:
    --------
        ses: list of sequences
        sequences for comparison

        ids: array of sequence id's

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
    opts, _ = getopt.getopt(argv, "a:m:s:g:", ['file', 'match', 'mismatch', 'gap'])
    
    sequences_array = []
    ids = []
    
    for opt, arg in opts:
        if opt in ["-a", "--file"]:

            #read in multiple sequences from fasta file
            input_file = SeqIO.parse(arg,'fasta')
            #sequences list
            sequences = [record for record in input_file]

            for record in sequences: 
                print (f"This is {record.description}: {record.seq}")
                #sequences_df = sequences.append(x)
                sequences_array.append(record.seq) 
                ids.append(record.description)

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
 
    return sequences_array,ids, match, mismatch, gap

def random_seq(L, rs):
    '''
    Parameters:
    ------------------
    L: (int) length of random sequence
    rs:(int) randomseed value

    Return:
    ----------------------
    random_dna_seq (String) random generated sequence
    counter (dic) value counts of bases in sequence
    '''
    random.seed(rs)
    dna = ["A","G","C","T"]
    random_dna_seq='' 
    
    for i in range(0,L):            
        random_dna_seq+=random.choice(dna)

    counter= Counter(random_dna_seq)
    return random_dna_seq, counter

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
            if str_seq0[row- 1] == "-" and str_seq1[column- 1] == "-": 
                vmatch = 0
            elif str_seq0[row- 1] == str_seq1[column- 1]:
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


def traceback_cross(S, T, rows, columns, sequence0, sequence1):
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
            tstring0 = tstring0 + "x"
            tstring1 = tstring1 + str_seq1[cur_col]
        else:
            # go one cell up
            cur_row = cur_row - 1 
            gap_no = gap_no + 1
            # add gap character to second sequence
            tstring0 = tstring0 + str_seq0[cur_row]
            tstring1 = tstring1 + "x"

    # reverse the completed alignment strings
    tstring0 = tstring0[::-1]
    tstring1 = tstring1[::-1]

    # get the optimal alignment score
    opt_score = S[rows-1, columns-1]
    return opt_score, match_no, mismatch_no, gap_no, tstring0, tstring1

    
def visual_alignment(tstring0, tstring1,filename, ids, match_no, mismatch_no, gap_no, match, mismatch, gap, pair_no=60 ):
    """
    Generate a visual alignment of two sequences following BLAST-alignment convention and write it to a text file
    (TASK 5)

    Parameters:
    -----------
        tstring0 (str): aligned string of first sequence
        tstring1 (str): aligned string of second sequence
        pair_no (int): number of aligned pairs to show in one row
        filename (String): tet.file name 
        id (): names of sequences
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
        file_out.write(f"\n ID of Sequence 1: {ids[0]} \n ID of Sequence 2: {ids[1]}")

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
      

def get_multiple_alignments(sequences,ids, match, mismatch, gap):
        
    # get matrices and number of rows/columns
    
    # get all alignment scores and alignment strings for combinations of sequences
    opt_scores = []
    list_strings0 = []
    list_strings1 = []
    for comb in list(map(list,combinations([0,1,2,3],2))):
        seqs_consider = itemgetter(comb[0], comb[1])(sequences)
        S, T, rows, columns =  compute(seqs_consider, match, mismatch, gap)
        opt_score, match_no, mismatch_no, gap_no, astring0, astring1 = traceback(S, T, rows, columns, seqs_consider[0], seqs_consider[1])
        # append combination strings and combination optimal scores
        list_strings0.append(astring0)
        list_strings1.append(astring1)
        opt_scores = np.append(opt_scores, opt_score)

    # get combinations into variable   
    combs = list(map(list,combinations([0,1,2,3],2)))
    # print("opt_scores = ", opt_scores)
    # print("List strings 0: ", list_strings0)

    # get the index of the optimal score
    max_score_arg = np.argmin(opt_scores)
    # get the aligned sequence strings with optimal scores
    A_max = [list_strings0[max_score_arg], list_strings1[max_score_arg]]
    print("A_max: ", A_max)
    #a_rest_inds = [x for x in [0,1,2,3] if x not in combs[max_score_arg]]
     
    # which sequences are NOT in A_max? Find out here
    rest_ind = [i for i, x in enumerate(combs) if x[0] not in combs[max_score_arg] and x[1] not in combs[max_score_arg]][0]
    # and put them into a rest
    A_rest = [list_strings0[rest_ind], list_strings1[rest_ind]]
    print("A_rest: ", A_rest)
    
    # in this part we get all possible cross-alignment-scores
    list_crossstrings0 = []
    list_crossstrings1 = []
    cross_opt_scores = []
    print(list(map(list,product([0,1],[0,1]))))
    for i, max_seq in enumerate(A_max):
        for j, rest_seq in enumerate(A_rest):
            S, T, rows, columns =  compute([max_seq, rest_seq], match, mismatch, gap)
            # in traceback_cross(), gaps are inserted as "x", so you can find the new gaps
            opt_score, match_no, mismatch_no, gap_no, astring0, astring1 = traceback_cross(S, T, rows, columns, max_seq, rest_seq)
            #print("Astring 0: ", astring0)
            cross_opt_scores = np.append(cross_opt_scores, opt_score)
            list_crossstrings0.append(astring0)
            list_crossstrings1.append(astring1)
            print("Optimal score for this combination: ", opt_score)

    # here you should find out the optimal cross score

    # and then you should insert gaps at the indices of the non-crossaligned strings where the cross-aligned strings have new gaps
    # this is not done yet
    max_cross_arg = np.argmin(cross_opt_scores)
    print("The cross combination with the best score: ", max_cross_arg)
    numpified_cross = np.array(list(str(list_crossstrings0[max_cross_arg])))
    print("Numpy array of strings of cross: ", np.array(list(str(list_crossstrings0[max_cross_arg]))))

    new_aligned = np.argwhere(numpified_cross == "x")

   
        
    print("Indices of new gaps: ", new_aligned)
    # combs_reduced = [x for x in combs if x not in [combs[max_score_arg], comb_rest_inds]]
    # inds_reduced = [i for i, x in enumerate(combs) if x not in [combs[max_score_arg], comb_rest_inds]]
    # opt_scores_reduced = [opt_scores[i] for i, x in enumerate(combs) if x not in [combs[max_score_arg], a_rest_inds]]

    
    #max_cross = np.argmax(opt_scores_reduced)
    #print("Max cross score: ", max_cross)
   # A_cross = [list_strings0[inds_reduced[max_cross]], list_strings1[inds_reduced[max_cross]]]
    #print("A_cross: ", A_cross)


    #visual_alignment(astring0, astring1, "dittschar_auckenthaler_assignment2_global_alignment.txt", ids, match_no, mismatch_no, gap_no, match, mismatch, gap)

def main():
     
    sequences,ids, match, mismatch, gap= get_args()
    # get matrices and number of rows/columns
    get_multiple_alignments(sequences,ids, match, mismatch, gap)

    #opt_scores = []
    # for comb in list(map(list,combinations([0,1,2,3],2))):
    #     seqs_consider = itemgetter(comb[0], comb[1])(sequences)
    #     S, T, rows, columns =  compute(seqs_consider, match, mismatch, gap)
    #     opt_score, match_no, mismatch_no, gap_no, astring0, astring1 = traceback(S, T, rows, columns, seqs_consider[0], seqs_consider[1])
    #     visual_alignment(astring0, astring1, "dittschar_auckenthaler_assignment2_global_alignment.txt", ids, match_no, mismatch_no, gap_no, match, mismatch, gap)
    #     opt_scores = np.append(opt_scores, opt_score)
        
    # combs = list(map(list,combinations([0,1,2,3],2)))
    # print("combinations: ", list(map(list,combinations([0,1,2,3],2))))
    # print("opt_scores = ", opt_scores)

    # max_score_arg = np.argmax(opt_scores)
    # A_max = [sequences[combs[max_score_arg][0]], sequences[combs[max_score_arg][1]]]
    
    # a_rest_inds = [x for x in [0,1,2,3] if x not in combs[max_score_arg]]
    # print("A rest inds: ", a_rest_inds)
    # A_rest = [sequences[a_rest_inds[0]], sequences[a_rest_inds[1]]]
    # combs_reduced = [x for x in combs if x not in [combs[max_score_arg], a_rest_inds]]
    # opt_scores_reduced = [opt_scores[i] for i, x in enumerate(combs) if x not in [combs[max_score_arg], a_rest_inds]]
    # print("Combs reduced: ", combs_reduced)
    # max_cross = np.argmax(opt_scores_reduced)
    # print("Max cross score: ", max_cross)
    #A_cross = 
    
    #S, T, rows, columns =  compute(sequences, match, mismatch, gap)
    # get optimal alignment score as well as number of matches, mismatches and gaps and aligned strings
    #opt_score, match_no, mismatch_no, gap_no, astring0, astring1 = traceback(S, T, rows, columns, sequences[0], sequences[1])
    #print("Optimal alignment score: ", opt_score)
    #call function to print and write output of needleman-wunsch
    #visual_alignment(astring0, astring1, "dittschar_auckenthaler_assignment2_global_alignment.txt", ids, match_no, mismatch_no, gap_no, match, mismatch, gap)

    



if __name__ == "__main__":
    try:
        main()
    except: 
        print("Try : python dittschar_auckenthaler_assignment2.py -a <file1> -m <match-score> -s <mismatch-score> -g <gap-score>")



    




    



