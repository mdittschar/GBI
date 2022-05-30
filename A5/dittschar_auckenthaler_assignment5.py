import getopt
import sys
import csv
import pandas as pd
import numpy as np
import itertools

def ccc(orig_dists, new_dists, num, no_files):
    """Computes CCC between two given dataframe matrices

    Args:
        orig_dists (pd.DataFrame): Original Distance Matrix
        new_dists (pd.DataFrame): Distance matrix of tree
        num        (int):           number of sequences
        no_files   (int):           number of files


    Returns:
        c (float): cophenetic correlation coefficient
    """
    orig_mean = np.mean(np.asarray(orig_dists).astype(float))
    new_mean = np.mean(np.asarray(new_dists).astype(float))
    numerator = []
    #number of files here relevant?
    numerator = np.append(numerator, 3* num)
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
    """get names and values of distance matrix

    Args:
        arg (file name): file name of matrix

    Returns:
        names (list of strings): sequence names
        values (array of floats): Values for each sequence
    """
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
                names.append(line[0])
                values.append(np.array(line_vals))       
    return names, values
        
def open_dist():
    """open dist files with given console arguments

    Returns:
        dist_df_1 (pd.DataFrame): first comparison dataframe
        dist_df_2 (pd.DataFrame): second comparison dataframe
        dist_df_3 (pd.DataFrame): third comparison dataframe
        mat_names (list of stirngs): File names
        num        (int):           number of sequences
        no_files   (int):           number of files

    """
    argv = sys.argv[1:]
    opts, _ = getopt.getopt(argv, "a:b:c:", ['file1', 'file2', 'file3'])
    
    names_list = []
    values_list = []
    mat_names = []
    for opt, arg in opts:
        if opt in ["-a", "--file1", "-b", "--file2", "-c", "--file3"]:
            print(f"File1: {arg}")
            names, values = get_mat(arg)
            names_list.append(names)
            values_list.append(values)
            mat_names.append(arg)
            
    
    dist_df_1 = pd.DataFrame(values_list[0], index=names_list[0])
    no_names_df1= len (dist_df_1.index)
   
    dist_df_2 = pd.DataFrame(values_list[1], index=names_list[1])
    no_names_df2= len (dist_df_2.index)
    dist_df_3 = pd.DataFrame(values_list[2], index=names_list[2])
    no_names_df3= len (dist_df_3.index)

    if (no_names_df1 & no_names_df2 == no_names_df3):
        num = no_names_df1
    no_files= (len(opts))

    return dist_df_1, dist_df_2, dist_df_3, mat_names, num, no_files



def main():
    orig_dist_df, dist_df_1, dist_df_2, mat_names, num, no_files = open_dist()
    c1 = ccc(orig_dists = orig_dist_df, new_dists=dist_df_1, num= num, no_files= no_files)
    c2 = ccc(orig_dists = orig_dist_df, new_dists=dist_df_2, num= num, no_files= no_files)
    print(f"CCC between {mat_names[0]} and {mat_names[1]} is: {c1}")
    print(f"CCC between {mat_names[0]} and {mat_names[2]} is: {c2}")




if __name__ == "__main__":
    try:
        main()
    except:
        print("Try python dittschar_auckenthaler_assignment5.py -a <original_distance> -b <first_tree_distance> -c <second_tree_distance>")