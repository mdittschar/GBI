import numpy as np
import sys
from Bio import SeqIO
import getopt


class Nussinov:
    def __init__(self):
        pass
    def matrix_filling(self, sequence, nmat):
        sequence = list(sequence)
        
        ptr = np.empty_like(nmat)
        ptr[:] = np.nan
        n = len(sequence)
        for m in np.arange(1, n):
            for j in np.arange(m, n):
                
                i = j - m #+ 1
                match = canonical(sequence[i], sequence[j])#[1 if sequence[i] == sequence[j] else 0][0]
                itrow = nmat[i+1, j]
                itcol = nmat[i, j-1]
                itdiag = nmat[i+1, j-1] + match
                bilist = [nmat[i, k] + nmat[k+1, j] for k in np.arange(i+1, j)]
                if len(bilist) > 0:
                    ititer = np.max(bilist)
                else:
                    ititer = 0


                nmat[i,j] = np.max(np.array([itrow, itcol, itdiag, ititer]))
        return nmat
        
    def matrix_init(self, sequence):
        sequence = list(sequence)
        n = len(sequence)
        nmat = np.empty((n, n))
        nmat[:] = np.nan
        for i in np.arange(n):
            nmat[i,i] = 0
            if i < n-1:
                nmat[i+1,i] = 0
        return nmat 
        
def traceback(nmat, sequence, structure, i, L):
    sequence = list(sequence)
    j = L
    if i < j:
        if nmat[i,j] == nmat[i+1, j]:
            traceback(nmat, sequence, structure, i+1, j)

        elif nmat[i,j] == nmat[i, j-1]:
            traceback(nmat, sequence, structure,i, j-1)
        elif nmat[i,j] == nmat[i+1, j-1] + canonical(sequence[i], sequence[j]):
            structure.append((i,j))
            traceback(nmat, sequence, structure,i+1, j-1)

        else:
            for k in np.arange(i+1, j-1):
                if nmat[i,j] == nmat[i,k] + nmat[k+1, j]:
                    traceback(nmat, sequence, structure, i, k)
                    traceback(nmat, sequence, structure, k+1, j)
                    break
    return structure
        
def dot_bracket(structure, sequence):
    basestring = ["." for i in np.arange(len(sequence))]
    
    for pair in structure:
        basestring[pair[0]] = "("
        basestring[pair[1]] = ")"
    print("Maximal number of basepairs: ", len(structure))
    return "".join(basestring)

def canonical(x,y):
    x = x[0]
    y = y[0]
    if x == "G" and y == "C":
        return 1
    elif x == "C" and y == "G":
        return 1
    elif x == "A" and y == "U":
        return 1
    elif x == "U" and y == "A":
        return 1
    else:
        return 0



def fasta_read():
    # fasta reader from console taken from our solution for A3
    argv = sys.argv[1:]
    # get a tuple of the arguments to use
    opts, _ = getopt.getopt(argv, "f:n:", ['file', 'name'])
    ses = [None] * 10

    i = 0
    for opt, arg in opts:
        if opt in ["-f", "--file"]:
            # use Biopython parser for sequences
            for record in SeqIO.parse(arg, "fasta"):
                se = record.seq
                # get sequences
                ses[i] = str(se)
                i = i + 1
    return ses



def main():
    # first sequence
    sequence = fasta_read()[0]
    print(f"Sequence: {sequence}")
    nuss = Nussinov()
    nmat = nuss.matrix_init(sequence)
    nmat = nuss.matrix_filling(sequence, nmat)
    structure = []
    structure = traceback(nmat, sequence, structure, 0, len(sequence)-1)
    dotbracketstring = dot_bracket(structure, sequence)
    print(dotbracketstring)
    with open(f"auckenthaler_dittschar_sequences.txt", "w") as f:   
        f.write(">sequence1")
        f.write("".join(sequence)+"\n")
        f.write(dotbracketstring)

    # second sequence
    sequence = fasta_read()[1]
    print(f"Sequence: {sequence}")
    nuss = Nussinov()
    nmat = nuss.matrix_init(sequence)
    nmat = nuss.matrix_filling(sequence, nmat)
    structure = list()
    sec = traceback(nmat, sequence, structure, 0, len(sequence)-1)
    print(structure)
    dotbracketstring = dot_bracket(structure, sequence)
    print(dotbracketstring)
    with open(f"auckenthaler_dittschar_sequences.txt", "a") as f:         
        f.write("\n>sequence2")
        f.write("".join(sequence)+"\n")
        f.write(dotbracketstring)

if __name__ == "__main__":
    main()
