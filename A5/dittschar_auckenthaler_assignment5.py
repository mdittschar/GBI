import getopt
import sys
import csv
import pandas as pd
import numpy as np

def ccc(orig_dists, new_dists):
    return

def get_df(arg):
        values = []
        names = []
        with open(arg) as file_in:
            for line in csv.reader(file_in, delimiter="\t"):    
                if len(line) > 1:                                  
                    names.append(line[0])
                    values.append(line[1:])
        dist_df = pd.DataFrame(values, index=names, columns=names)
def open_dist():
    argv = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "a:b:c:", ['file1', 'file2', 'file3'])
    for opt, arg in opts:
        if opt in ["-a", "--file1", "-b", "--file2", "-c", "--file3"]:
            print(f"File1: {arg}")
            dist_df = get_df(arg)

            
    return dist_df



def main():
    open_dist()
    #ccc(orig_dists, new_dists)




if __name__ == "__main__":
    try:
        main()
    except:
        print("Try python dittschar_auckenthaler_assignment5.py")