import argparse
import numpy as np
import math
from itertools import combinations_with_replacement


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

# aus dem Internet kopiert, noch ändern!
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

def fasta_len(seq):
    print(len(seq))
    return len(seq)

def rev_complement(seq):
    
    rev_seq = seq[::-1]
    list_seq = list(rev_seq)
    rev_comp = []
    for char in list_seq:
        if char == 'A':
            rev_comp = np.append(rev_comp, 'T')
        elif char == 'G':
            rev_comp = np.append(rev_comp, 'C')
        elif char == 'T':
            rev_comp = np.append(rev_comp, 'A')
        elif char == 'C':
            rev_comp = np.append(rev_comp, 'G')

    return ''.join(rev_comp)

def read_file(msf):
    names = []
    seqs_list = []
    with open(msf) as fp:
        for name, seq in read_fasta(fp):
            seq = list(seq)
            seqs_list = seqs_list +  seq
            names = np.append(names, name)
    return names, seqs_list

def get_combinations(seqs2, names2):
    n=int(len(seqs2)/len(names2))

    seqs_matrix=[seqs2[i:i + n] for i in range(0, len(seqs2), n)]
    seqs_matrix = np.array(seqs_matrix)
    print("Output list chunks: ", seqs_matrix)

    # seqs_str = ""

    #print(seqs_str)
    # with open(file_out, 'w') as f_out:   
    #     f_out.write(">" + name + "\n" + seq + "\n")         
        
   

        
    
    # T3

    # pass is used when a function is not completely difened.
    # Remove it once you have started with the assignnment.

    

    # s1= "AAABCAABBC"
    # s2= "ABABCBABBB"
    # s3= "AACBCCABBA"
    # #convert to list
    # seq1= list(s1)
    # seq2= list(s2)
    # seq3= list(s3)
    
    # # #Add to one sequence
    # seqs= seq1+seq2+seq3
    
    
    # #print(len(seq3))
    # print("Type of seqs: ", type(seqs))

    unique = dict(zip(seqs2,[seqs2.count(i) for i in seqs2]))
    #rint("Unique with string: ", unique)

    unique_val= np.unique(np.array(list(seqs2)))

    combis= (list(combinations_with_replacement(unique_val, 2))) 
    combis= np.array(combis)
    # seqs_matrix= [seq1, seq2, seq3]
    
    seqs_matrix= np.array(seqs_matrix)
    #print("Sequences matrix", seqs_matrix)

    combi_counts=dict()
    for combi in range (len(combis)):
        count = 0
        for column in np.arange(seqs_matrix.shape[1]):
            for row in np.arange(seqs_matrix.shape[0]-1):
                for compare in np.arange(row, seqs_matrix.shape[0]):
                    if compare != row: 
                        c = ()
                        c0= seqs_matrix[row,column]
                        c1= seqs_matrix[compare,column]
                        #print("C0: ", c0, "C1: ", c1)
                        cc= [[c0,c1]]
                        #print(c0,c1)
                        if (combis[combi,0]!= combis[combi,1]):
                            if (((c0 == combis[combi,0]) & (c1 == combis[combi,1])) | ((c0 == combis[combi,1]) & (c1 == combis[combi,0]))):
                                count = count + 1
                        else:
                            if ((c0 == combis[combi,0] ) & (c1 == combis[combi,1])):
                                count = count+ 1   
                         
                    else: 
                        pass
            #print("Looking for", combis[combi,:], "counts: ", count)
        combi_counts[(np.array_str(combis[combi,:]))] = count
    print(combi_counts)
    return combi_counts, unique_val


def main():
    '''
    The main function should contain all functions that solved the stated tasks. 
    '''
    # T2.a
    msf = args.file_one
    file_out='gene_seq_out.fasta'
    msf2 = args.file_two
    names, seqs = read_file(msf)  
    names2, seqs2 = read_file(msf2)     

    # print("No of names" ,len(names2))     
    # print("Lenght of seqs: " ,len(seqs2))

    # n=int(len(seqs2)/len(names2))

    # seqs_matrix=[seqs2[i:i + n] for i in range(0, len(seqs2), n)]
    # seqs_matrix = np.array(seqs_matrix)
    

    
    combi_counts, unique_val = get_combinations(seqs2, names2)

  

    # S =np.empty((len(unique_val),len(unique_val)))
    # combi_counts_keys= list(combi_counts.keys())
    # unique_counts_keys = list(unique.keys())
    # S[0,0] = math.log2((combi_counts[combi_counts_keys[0]]/30)/((unique[unique_counts_keys[0]]/len(seqs))*(unique[unique_counts_keys[0]]/len(seqs))))
    # S[0,1] =S[1,0]= math.log2((combi_counts[combi_counts_keys[1]]/30)/((unique[unique_counts_keys[0]]/len(seqs))*(unique[unique_counts_keys[1]]/len(seqs))))
    # S[1,1] = math.log2((combi_counts[combi_counts_keys[3]]/30)/((unique[unique_counts_keys[1]]/len(seqs))*(unique[unique_counts_keys[1]]/len(seqs))))
    # S[0,2] =S[2,0] = math.log2((combi_counts[combi_counts_keys[2]]/30)/((unique[unique_counts_keys[0]]/len(seqs))*(unique[unique_counts_keys[2]]/len(seqs))))
    # S[2,2] = math.log2((combi_counts[combi_counts_keys[5]]/30)/((unique[unique_counts_keys[2]]/len(seqs))*(unique[unique_counts_keys[2]]/len(seqs))))
    # S[1,2] =S[2,1] = math.log2((combi_counts[combi_counts_keys[4]]/30)/((unique[unique_counts_keys[1]]/len(seqs))*(unique[unique_counts_keys[2]]/len(seqs))))
    # print("S: ", S)

    pass


if __name__ == "__main__":
    try:
        args = create_parser()
        # accesing the path of the files
        print(args.file_one)
        print(args.file_two)

        main()
    except:
        print('Try:  python template-a1.py -f1 MultipleSeqs.fasta -f2 msa-scoring-matrix.fasta')

