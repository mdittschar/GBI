import getopt
import sys
import numpy as np

def get_seqs():
    """From a fasta file, get sequences and ids
    Taken from our previous assignments
    Returns:
        astring (list): Letters from the first input (PHD)
        bstring (list): Letters from the second input (Proteus2)
        rstring (list): Letters from the third input (reference)
    """
    argv = sys.argv[1:]
    # get a tuple of the arguments to use
    opts, _ = getopt.getopt(argv, "a:b:r:", ['file1','file2', 'reference'])
    for opt, arg in opts:
        if opt in ["-a", "--file1"]:        
            with open(arg) as f:
                astring = f.readlines()
                astring = list(astring[0])
        elif opt in ["-b", "--file2"]:        
            with open(arg) as f:
                bstring = ""
                for i, line in enumerate(f.readlines()):
                    if i %5 == 2: # parsing is adapted to input file configuration
                        bstring = bstring + line
                bstring = list(bstring)
        elif opt in ["-r", "--reference"]:        
            with open(arg) as f:
                rstring = ""
                for i, line in enumerate(f.readlines()):
                    if i %4 == 2: # parsing is adapted to input file configuration
                        rstring = rstring + line[8:]
                rstring = list(rstring)

    return astring, bstring, rstring

def equalize_strings(astring):
    """Transform letters in a list into broader structure categories

    Args:
        astring (list): Letters

    Returns:
        astring (list): transformed letters
    """
    astring = ["H" if i in ["H", "G", "I"] else i for i in astring]
    astring = ["E" if i in ["E", "B"] else i for i in astring ]
    astring = ["C" if i in ["T", "S", "C", "L"] else i for i in astring ]
    return astring

def compare_strings(astring, bstring, rstring):
    """Return Q3-scoretwo different strings  and a reference

    Args:
        astring (list): phd structure letters
        bstring (list): proteus2 structure letters
        rstring (list): reference structure letters
    """
    # remove all symbols that aren't true letters
    astring = [i for i in astring if i in ["E", "H", "L", "C"]]
    bstring = [i for i in bstring if i in ["E", "H", "L", "C"]]
    rstring = [i for i in rstring if i in ["E", "H", "C", "T", "G", "B"]]

    astring = equalize_strings(astring)
    bstring = equalize_strings(bstring)
    rstring = equalize_strings(rstring)
    astring = np.array(astring)
    bstring = np.array(bstring)
    rstring = np.array(rstring)

    numerator_a = np.sum(astring==rstring)
    numerator_b = np.sum(bstring==rstring)
    print(f"Q3 of PHD sequence is: \n{numerator_a/len(astring)}")
    print(f"Q3 of Proteus 2 sequence is: \n{numerator_b/len(bstring)}")

def main():
    astring, bstring, rstring = get_seqs()
    compare_strings(astring, bstring, rstring)


if __name__ == "__main__":
    main()