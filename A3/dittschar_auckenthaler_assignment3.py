import Bio
from Bio import SeqIO
import getopt
import sys
import numpy as np
import pandas as pd
from itertools import combinations, product
from operator import itemgetter
import math
import random
from collections import Counter

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
    opts, _ = getopt.getopt(argv, "a:b:cdm:s:g:", ['file1', 'file2', 'file3', 'file4', 'match', 'mismatch', 'gap'])
    ses = [None] * 4

    ids = []
    i = 0
    for opt, arg in opts:
        if opt in ["-a", "--file1", "-b", "--file2", "-c", "--file3", "-d", "--file4"]:
            print(f"File No {i+1}: {arg}")
            # use Biopython parser for sequences
            for record in SeqIO.parse(arg, "fasta"):
                print("record: ", i)
                se = record.seq
                
                id = record.id
                
                #ids [i]= id
                # append sequences to list
                ids = np.append(ids, id)

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
 
    return ses,ids, match, mismatch, gap

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

    
    

def traceback(S, T, rows, columns, sequence0, sequence1 ,gap_c):
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
        gap_c (String):   gap character

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
            tstring0 = tstring0 + gap_c
            tstring1 = tstring1 + str_seq1[cur_col]
        else:
            # go one cell up
            cur_row = cur_row - 1 
            gap_no = gap_no + 1
            # add gap character to second sequence
            tstring0 = tstring0 + str_seq0[cur_row]
            tstring1 = tstring1 + gap_c

    # reverse the completed alignment strings
    tstring0 = tstring0[::-1]
    tstring1 = tstring1[::-1]

    # get the optimal alignment score
    opt_score = S[rows-1, columns-1]
    return opt_score, match_no, mismatch_no, gap_no, tstring0, tstring1

  
def visual_alignment(tstring0, tstring1, tstring2, tstring3, filename, ids, match, mismatch, gap, pair_no=60 ):
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
        #for k in range(len(ids)):
            #file_out.write(f"\n ID of Sequence {k+1}: {ids[k]}")

        while i*pair_no < seq_lens:
            # assign strings to be visualised
            # pair_no is number of pairs to visualise in one row
            display_seq0 = tstring0[i*pair_no:(i+1)*pair_no]
            display_seq1 = tstring1[i*pair_no:(i+1)*pair_no]
            display_seq2 = tstring2[i*pair_no:(i+1)*pair_no]
            display_seq3 = tstring3[i*pair_no:(i+1)*pair_no]
            # show alignment lines where matches are present
            match_array0 = np.where(np.logical_and(np.array(list(tstring0))== np.array(list(tstring1)),np.array(list(tstring0)) != "-"), "|", " ")
            match_array1 = np.where(np.logical_and(np.array(list(tstring1))== np.array(list(tstring2)),np.array(list(tstring1)) != "-"), "|", " ")
            match_array2 = np.where(np.logical_and(np.array(list(tstring2))== np.array(list(tstring3)),np.array(list(tstring2)) != "-"), "|", " ")
            no_gaps_array0 = np.logical_and(np.array(list(tstring0)) != "-", np.array(list(tstring1)) != "-")
            no_gaps_array1 = np.logical_and(np.array(list(tstring1)) != "-", np.array(list(tstring2)) != "-")
            no_gaps_array2 = np.logical_and(np.array(list(tstring2)) != "-", np.array(list(tstring3)) != "-")
            display_alignment0 = "".join(np.where(np.logical_and(no_gaps_array0,np.array(list(tstring0)) != np.array(list(tstring1))), "*", match_array0)[i*pair_no:(i+1)*pair_no])
            display_alignment1 = "".join(np.where(np.logical_and(no_gaps_array1, np.array(list(tstring1)) != np.array(list(tstring2))), "*", match_array1)[i*pair_no:(i+1)*pair_no])
            display_alignment2 = "".join(np.where(np.logical_and(no_gaps_array2,np.array(list(tstring2)) != np.array(list(tstring3))), "*", match_array2)[i*pair_no:(i+1)*pair_no])

            #display_alignment0 = "".join(np.where(np.array(list(tstring0)) == '-', " ", match_array0)[i*pair_no:(i+1)*pair_no])
            #display_alignment1 = "".join(np.where(np.array(list(tstring1)) == '-', "i ", match_array1)[i*pair_no:(i+1)*pair_no])
            #display_alignment2 = "".join(np.where(np.array(list(tstring2)) == '-', " i", match_array2)[i*pair_no:(i+1)*pair_no])

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
            print(f"\n{alignment_nos}\n{display_seq0}\n{display_alignment0}\n{display_seq1}\n{display_alignment1}\n{display_seq2}\n{display_alignment2}\n{display_seq3}")
            #Task 5 a) write it to textfile
            file_out.write(f"\n{alignment_nos}\n{display_seq0}\n{display_alignment0}\n{display_seq1}\n{display_alignment1}\n{display_seq2}\n{display_alignment2}\n{display_seq3}")
        
            i = i + 1
        #Task 5 b) and c) write parameter to text file:
        print("Match Score: ",match)
        print("Mismatch Score: ", mismatch)
        print("Gap penalty: ",gap )
       
        file_out.write(f"\nMatch Score: {match}\nMismatch Score: {mismatch}\nGap Penalty: {gap}")
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
        opt_score, match_no, mismatch_no, gap_no, astring0, astring1 = traceback(S, T, rows, columns, seqs_consider[0], seqs_consider[1], '-')
        # append combination strings and combination optimal scores
        list_strings0.append(astring0)
        list_strings1.append(astring1)
        opt_scores = np.append(opt_scores, opt_score)

    
    # get combinations into variable   
    combs = list(map(list,combinations([0,1,2,3],2)))
    # get the index of the optimal score
    max_score_arg = np.argmax(opt_scores)
    # get the aligned sequence strings with optimal scores
    A_max = [list_strings0[max_score_arg], list_strings1[max_score_arg]]
    #a_rest_inds = [x for x in [0,1,2,3] if x not in combs[max_score_arg]]
     
    # which sequences are NOT in A_max? Find out here
    rest_ind = [i for i, x in enumerate(combs) if x[0] not in combs[max_score_arg] and x[1] not in combs[max_score_arg]][0]
    # and put them into a rest
    A_rest = [list_strings0[rest_ind], list_strings1[rest_ind]]
    
    
    # in this part we get all possible cross-alignment-scores
    list_crossstrings0 = []
    list_crossstrings1 = []
    cross_opt_scores = []

    combis_Amax_Arest= (list(map(list,product([0,1],[0,1]))))
    for i, max_seq in enumerate(A_max):
        for j, rest_seq in enumerate(A_rest):
            S, T, rows, columns =  compute([max_seq, rest_seq], match, mismatch, gap)
            opt_score, match_no, mismatch_no, gap_no, astring0, astring1 = traceback(S, T, rows, columns, max_seq, rest_seq,'x')
            cross_opt_scores = np.append(cross_opt_scores, opt_score)
            list_crossstrings0.append(astring0)
            list_crossstrings1.append(astring1)
            

    #find out the optimal cross score 
    max_cross_arg = np.argmax(cross_opt_scores)

    A_max_cross = np.array(list(str(list_crossstrings0[max_cross_arg])))

    A_rest_cross = np.array(list(str(list_crossstrings1[max_cross_arg])))

    idx_of_x_max = np.ndarray.flatten(np.argwhere(A_max_cross == "x"))
    idx_of_x_rest = np.ndarray.flatten(np.argwhere(A_rest_cross == "x"))

    
    if combis_Amax_Arest[max_cross_arg][0] == 1:
        A_max_noncross = A_max[0]
    else: 
        A_max_noncross = A_max[1]

    if combis_Amax_Arest[max_cross_arg][1] == 1:
        A_rest_noncross = A_rest[0]
    else: 
        A_rest_noncross = A_rest[1]

    for i in idx_of_x_max:
        A_max_noncross = list(A_max_noncross)
        A_max_noncross.insert(i, "x")
    for i in idx_of_x_rest:
        A_rest_noncross = list(A_rest_noncross)
        A_rest_noncross.insert(i, "x")
        

    A_rest_noncross = np.array(A_rest_noncross)
    A_max_noncross = np.array(A_max_noncross)
    A_rest_cross = np.array(A_rest_cross)
    A_max_cross = np.array(A_max_cross)

    A_max_cross = np.where(A_max_cross =="x", "-", A_max_cross)
    A_max_noncross = np.where(A_max_noncross =="x", "-", A_max_noncross)
    A_rest_cross = np.where(A_rest_cross =="x", "-", A_rest_cross)
    A_rest_noncross = np.where(A_rest_noncross =="x", "-", A_rest_noncross)

    A_max_cross = "".join(A_max_cross)
    A_max_noncross = "".join(A_max_noncross)
    A_rest_cross = "".join(A_rest_cross)
    A_rest_noncross = "".join(A_rest_noncross)


    return  A_max_cross, A_max_noncross, A_rest_cross, A_rest_noncross
    


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

def feng_doolittle_distance(sequence0, sequence1, match, mismatch, gap, L):
    '''
    Parameters:
    ----------------
    sequence0:(String)  sequence 1 
    sequence1:(String)  seqeuence 2
    match:(int)         match score
    mismatch:(int)      mismatch score
    gap: (int)          gap score
    L: (Int)            Length of sequences

    Return:
    --------------
    d = (float)         calculated feng-doolittle-distance
    '''
    #S_obs (X,Y)
    S_obs_seqs= [sequence0,sequence1]
    S_obs, T_obs, rows_obs, columns_obs =  compute(S_obs_seqs, match, mismatch, gap)
    S_obs, match_no, mismatch_no, gap_no, astring0, astring1 = traceback(S_obs, T_obs, rows_obs, columns_obs, S_obs_seqs[0], S_obs_seqs[1],'-')
    #print("Feng- Doolittle S_obs: ", S_obs)
    #Optimal Score S(X,X)
    S_id0_seqs= [sequence0,sequence0]
    S_id0, T_id0, rows_id0, columns_id0 =  compute(S_id0_seqs, match, mismatch, gap)
    S_id0_ops, match_no__id0, mismatch_no__id0, gap_no__id0, astring0__id0, astring1__id0 = traceback(S_id0, T_id0, rows_id0, columns_id0, S_id0_seqs[0], S_id0_seqs[1],'-')
    #optimal Score S(Y,Y)
    S_id1_seqs= [sequence1,sequence1]
    S_id1, T_id1, rows_id1, columns_id1 =  compute(S_id1_seqs, match, mismatch, gap)
    S_id1_ops, match_no__id1, mismatch_no__id1, gap_no__id1, astring0__id1, astring1__id1 = traceback(S_id1, T_id1, rows_id1, columns_id1, S_id1_seqs[0], S_id1_seqs[1],'-')
    #S_id
    S_id= (S_id0_ops+S_id1_ops)/2
    #print("Feng- Doolittle S_id: ",S_id)

    random_seqX, counterX= random_seq(L, 1)
    random_seqY, counterY= random_seq(L, 2)

    Ns= []
    for a in (["A","G","C","T"]):
        NX= counterX[a]
        #print (NX)
        NY= counterY[a]
        #print(NY)
        Ns = np.append(Ns,(NX*NY*mismatch))
        #print("N von a: ", Ns)
    #print(Ns)
    N= sum(Ns)
    #print (N)
    S_rand_seqs= [random_seqX, random_seqY]
    #print(S_rand_seqs)
    S_rand_S, T_rand, rows_rand, columns_rand =  compute(S_rand_seqs, match, mismatch, gap)
    S_rand_obs, match_no_rand, mismatch_no_rand, gap_no_rand, astring0, astring1 = traceback(S_rand_S, T_rand, rows_rand, columns_rand, S_rand_seqs[0], S_rand_seqs[1],'-')

    S_rand= 1/L *N - gap_no_rand*gap
    #print("Feng- Doolittle S_rand: ",S_rand)
    d= -math.log((S_obs-S_rand)/(S_id- S_rand))
    #print("Distance: ",d)
    return d

def distance_matrix(sequences,ids, filename, match, mismatch, gap, L): 
    '''
    Parameters:
    ------------
    sequences [array]   array of sequences from fasta file
    filename [String]   name of the file to store distance matrix
    match (int)         match-score
    mismatch (int)      mismatch score
    gap (int)           gap-score
    L(int)              lenth of random sequence 

    Return:
    --------------
    d_matrix_df(Dataframe)  Distance matrix    
    '''
    with open(filename, 'w') as file_out:
        d_matrix = np.zeros((len(sequences), len(sequences)))
        for j in range(len(sequences)):
            for i in range(len(sequences)):
                seqX= sequences[i]
                seqY= sequences[j]
                d= feng_doolittle_distance(seqX,seqY,match,mismatch, gap,L)

                d_matrix[i][j]= d
        d_matrix_df= pd.DataFrame(d_matrix, index=[1,2,3,4], columns=[1,2,3,4])
        print (f"Distance Matrix: \n{d_matrix_df}")
        #for k in range(len(ids)):
            #file_out.write(f"\n ID of Sequence {k+1}: {ids[k]}")
        file_out.write(f"\nDistance Matrix:\n {d_matrix_df}")
        file_out.close


    return d_matrix_df

def main():
     
    sequences,ids, match, mismatch, gap= get_args()
    # get matrices and number of rows/columns
    profile1_0, profile1_1, profile2_0, profile2_1 = get_multiple_alignments(sequences,ids, match, mismatch, gap)

    

    distance_matrix(sequences,ids,"dittschar_auckenthaler_assignment3_distance_matrix.txt",match, mismatch, gap, L=60)
    
    #d= feng_doolittle_distance(sequence0, sequence1, match, mismatch, gap, L= 60)
      
   
    
    
    #S, T, rows, columns =  compute(sequences, match, mismatch, gap)
    # get optimal alignment score as well as number of matches, mismatches and gaps and aligned strings
    #opt_score, match_no, mismatch_no, gap_no, astring0, astring1 = traceback(S, T, rows, columns, sequences[0], sequences[1])
    #print("Optimal alignment score: ", opt_score)
    #call function to print and write output of needleman-wunsch
    visual_alignment(profile1_0, profile1_1, profile2_0, profile2_1, "dittschar_auckenthaler_assignment3_profile_alignment.txt", ids, match, mismatch, gap)

    



if __name__ == "__main__":
    try:
        main()
    except: 
        print("Try : python dittschar_auckenthaler_assignment2.py -a <file1> -b <file2> -m <match-score> -s <mismatch-score> -g <gap-score>")

