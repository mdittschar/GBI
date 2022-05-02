import argparse
import numpy as np


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
    seqs = []
    with open(msf) as fp:
        for name, seq in read_fasta(fp):
            seqs = np.append(seqs, seq)
            names = np.append(names, name)
    return names, seqs


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

    print(names)
    print(seqs.shape)      
    print(np.expand_dims(seqs, axis=0))
    # Man kann einzelne Buchstaben ausgeben lassen!
    print(np.expand_dims(seqs, axis=0)[-1][-1])

    #print(seqs[0,0] == 'A')
    #print(rev_complement(seqs[-1]))
    for seq in np.asarray(list(seqs)):
        list_seq = np.asarray(list(seq))
        print(np.sum(list_seq=='A'))
    # with open(file_out, 'w') as f_out:   
    #     f_out.write(">" + name + "\n" + seq + "\n")         
        
        

        
    
    # T3

    # pass is used when a function is not completely difened.
    # Remove it once you have started with the assignnment.
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

