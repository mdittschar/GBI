import argparse


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


def main():
    '''
    The main function should contain all functions that solved the stated tasks. 
    '''
    # T2.a

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
        print('Try:  python3 template-a1.py -f1 MultipleSeqs.fasta -f2 msa-scoring-matrix.fasta')
