import argparse
import numpy as np
import math
from math import comb, log2
from itertools import combinations_with_replacement
import pandas as pd
from collections import OrderedDict

def create_parser():
    '''
    Creates a argument parser to facilitate file reading.

    '''
    # Make parser object
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('-f1', '--file-one',
                   help="File for the the first programming task")
    p.add_argument('-f2', '--file-two',
                   help="File for the the second programming task")

    return(p.parse_args())


def get_name_seq(f):
    '''
    Seperate names 

    Parameters
    ----------
    f:file , (opened file)
    '''
    name = False
    seq = []
    for line in f:
        line = line.rstrip()
        if line[0] == ">":
            if name: 
                yield (name, ''.join(seq))
            name= line
            seq = []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def fasta_len(seq):
    '''
    Calculates the length of a sequence

    Parameters
    ----------
    seq: list or np.array, (sequence)

    Returns
    ----------
    length: int
    length of sequence
    '''
    length= len(seq)
    print(length)
    return length

def rev_complement(seq):
    '''
    Generate reverse compliment of a sequence

    Parameters
    ----------
    seq: list or np.arry, (sequence)
    Sequence

    Returns
    ----------
    rev_comps: list, (n_reversed_compliment_bases)
    Reverse compliment sequence
    '''
    # reverse the sequence
    rev_seq = seq[::-1]
    list_seq = list(rev_seq)
    rev_comp = []
    # replace by complement
    for char in list_seq:
        if char == 'A':
            rev_comp = np.append(rev_comp, 'T')
        elif char == 'G':
            rev_comp = np.append(rev_comp, 'C')
        elif char == 'T':
            rev_comp = np.append(rev_comp, 'A')
        elif char == 'C':
            rev_comp = np.append(rev_comp, 'G')
    rev_comp_seq= ''.join(rev_comp)
    return rev_comp_seq

def read_file(msf):
    '''
    read sequences from file

    Parameters
    ----------
    msf: file, (fasta file)

    Returns
    ----------
    names: np.array, (n_sequence names)
    seqs: np.array, (n_sequences)
    seqs_list: list, (sequence)
    '''

    try:
        names = []
        seqs_list = []
        seqs = []
        with open(msf) as f:
            for name, seq in get_name_seq(f):
                seq_list = list(seq)
                seqs_list = seqs_list +  seq_list
                names = np.append(names, name)
                seqs = np.append(seqs, seq)
        return names, seqs, seqs_list
    except FileNotFoundError:
        print("File not found!")


def write_file(path, names, seqs):
    '''
    write a file

    Parameters
    ----------
    path: string, 
    names: np.array (n_names,)
    seqs: np.array (n_sequences)
    '''
    with open(path, 'w') as f_out:
        for name, seq in zip(names, seqs):

            f_out.write(">" + name + "\n" + seq + "\n")


def get_combinations(seqs2, names2):
    '''
    Get and count all possible combinations, and single values in sequence

    Parameters
    ----------
    seqs2: list of sequences
    sequences
    names2: np.array
    names

    Returns
    ----------
    combi_counts: dic, (combi, count) 
    Counts for every combbination
    unique: dic, (c, count)
    probability of single letters

    '''
    # get the length of each sequence
    n=int(len(seqs2)/len(names2))

    # generate a matrix from big list
    seqs_matrix=[seqs2[i:i + n] for i in range(0, len(seqs2), n)]
    seqs_matrix = np.array(seqs_matrix)
     
    
    # get unique counts
    unique = OrderedDict(sorted(dict(zip(seqs2,[seqs2.count(i) for i in seqs2])).items()))
    
    # get unique letters
    unique_val= np.unique(np.array(list(seqs2)))
    # generate all possible combinations
    combis= (list(combinations_with_replacement(unique_val, 2))) 
    combis= np.array(combis)
    
    seqs_matrix= np.array(seqs_matrix)

    # for each combination, for each column and row and comparison row, get combinations
    combi_counts=dict()
    for combi in range (len(combis)):
        count = 0
        for column in np.arange(seqs_matrix.shape[1]):
            for row in np.arange(seqs_matrix.shape[0]-1):
                for compare in np.arange(row, seqs_matrix.shape[0]):
                    if compare != row: 
                        c0= seqs_matrix[row,column]
                        c1= seqs_matrix[compare,column]
                        if (combis[combi,0]!= combis[combi,1]):
                            # conditions for letter comparison
                            if (((c0 == combis[combi,0]) & (c1 == combis[combi,1])) | ((c0 == combis[combi,1]) & (c1 == combis[combi,0]))):
                                count = count + 1
                        else:
                            if ((c0 == combis[combi,0] ) & (c1 == combis[combi,1])):
                                count = count+ 1   
                         
                    else: 
                        pass
        combi_counts[(np.array_str(combis[combi,:]))] = count
    
    return combi_counts, unique

def compute_sub_matrix(combi_counts, unique, seqs2_list):
    '''
    Calculate Substitution Matrice

    Parameters
    ----------
    combi_counts: dic, (combi, count) 
    unique: dic, (c, count)
    seqs2_list: list , (sequences)

    Returns
    ----------
    df: pd.Dataframe, (nxn_matrix)

    '''
    seq_len = len(seqs2_list)
    # generate dataframe
    df = pd.DataFrame(unique, index=unique.keys())
    no_all_combis = sum(combi_counts.values())
    # get two letters
    for key in combi_counts.keys():
        # the letters are always at the same position
        count1 = key[2:3]
        count2 = key[6:7]
        # p(A,B)/(p(A)*p(B))
        val = (combi_counts[key]/no_all_combis)/((unique[count1]/seq_len)*(unique[count2]/seq_len))
        # put in values
        if val ==0:
            df.loc[count1, count2] = 0
            df.loc[count2, count1] = 0
        else:
            logval = math.log(val,2)
            df.loc[count1, count2] = np.round(logval, 2)
            df.loc[count2, count1] = np.round(logval, 2)
        
    return df

def matrix_to_txt(matrix, path):
    """
    write matrix to text file

    Parameters
    ----------

    matrix: pd.DataFrame
    (Substitution) Matrix

    path: string 
    filepath
    """
    matrix.to_csv(path, sep='\t', mode='a')

def main():
    '''
    The main function should contain all functions that solved the stated tasks. 
    '''
    # Task 2
    msf = args.file_one
    file_out='seq_out.fasta'
    msf2 = args.file_two
    names, seqs, seqs_list = read_file(msf)  
    names2, seqs2, seqs2_list= read_file(msf2)     
    write_file(path=file_out, names=names, seqs=seqs)
    # Task 3
    combi_counts, unique = get_combinations(seqs2_list, names2)
    matrix = compute_sub_matrix(combi_counts, unique, seqs2_list)
    print("Substitution matrix: \n", matrix)
    matrix_to_txt(matrix, path="matrix.txt")
    


if __name__ == "__main__":
    try:
        args = create_parser()
        # accesing the path of the files
        print(args.file_one)
        print(args.file_two)

        main()
    except:
        print('Try:  python template-a1.py -f1 MultipleSeqs.fasta -f2 msa-scoring-matrix.fasta')

