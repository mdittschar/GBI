import getopt
import sys
import csv
import pandas as pd
import numpy as np
import itertools

def ccc(orig_dists, new_dists):
    orig_mean = np.mean(np.asarray(orig_dists).astype(float))
    new_mean = np.mean(np.asarray(new_dists).astype(float))
    numerator = []
    numerator = np.append(numerator, 3*9)
    denominator_one = []
    denominator_two = []
    for i, val in enumerate(orig_dists.index):
        for j in np.arange(i+1):
            numerator_one = orig_dists.loc[val, j] - orig_mean
            numerator_two = new_dists.loc[val, j] - new_mean
            
            numerator = np.append(numerator, numerator_one * numerator_two)
            denominator_one = np.append(denominator_one, numerator_one * numerator_one)
            denominator_two = np.append(denominator_two, numerator_two * numerator_two)
    numerator_sum = sum(numerator)
    denominator_sum_one = sum(denominator_one)
    denominator_sum_two = sum(denominator_two)

    c = numerator_sum/(np.sqrt(denominator_sum_one*denominator_sum_two))
    return c

def get_mat(arg):
        values = []
        names = []
        with open(arg) as file_in:
            for idx, line in enumerate(csv.reader(file_in, delimiter="\t")):
                if idx >0:
                    if idx == 1:
                        first_val = line[1]
                        line_vals = [float(first_val)] + [0]*(9-idx)
                    else: 
                        line_vals = [float(i) for i in line[1:]] + [0]*(9-idx)
                    print("Name: ",line[0])                             
                    names.append(line[0])
                    
                    print(line_vals)
                    values.append(np.array(line_vals))
        return names, values
        
def open_dist():
    argv = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "a:b:c:", ['file1', 'file2', 'file3'])
    names_list = []
    values_list = []
    for opt, arg in opts:
        if opt in ["-a", "--file1", "-b", "--file2", "-c", "--file3"]:
            print(f"File1: {arg}")
            names, values = get_mat(arg)
            names_list.append(names)
            values_list.append(values)
            
    
    dist_df_1 = pd.DataFrame(values_list[0], index=names_list[0])
    dist_df_2 = pd.DataFrame(values_list[1], index=names_list[1])
    dist_df_3 = pd.DataFrame(values_list[2], index=names_list[2])
    
    return dist_df_1, dist_df_2, dist_df_3



def main():
    orig_dist_df, dist_df_1, dist_df_2 = open_dist()
    print(dist_df_2)
    c1 = ccc(orig_dists = orig_dist_df, new_dists=dist_df_1)
    c2 = ccc(orig_dists = orig_dist_df, new_dists=dist_df_2)
    print(f"C1 is: {c1}\nC2 is {c2}")




if __name__ == "__main__":
    try:
        main()
    except:
        print("Try python dittschar_auckenthaler_assignment5.py")