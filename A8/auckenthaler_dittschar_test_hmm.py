import numpy as np
import sys
import getopt
from Bio import SeqIO
import pandas as pd
from io import StringIO
from collections import OrderedDict

argv = sys.argv[1:]
# get a tuple of the arguments to use
opts, _ = getopt.getopt(argv, "a:b:i:", ['matrix1', 'matrix2', 'input'])


for opt, arg in opts:
    if opt in ["-a", "--matrix1"]:
        with open(arg) as f:
            lines = f.readlines()
        # get the states
        states = lines[3]
        # every second element in string is a state
        states = list(states)[::2]
        print("States:", states)
        # exclude the lines that are not transition matrix
        lines = lines[-4:]
        plusmat = np.empty((len(lines), len(lines)))
        # convert from a text matrix to a numpy array
        for i, line in enumerate(lines): 
            plusmat[i, :] = np.squeeze(np.array(pd.read_csv(StringIO(line[:-1]),header=None,  sep=" "))).astype(float)

    if opt in ["-b", "--matrix2"]:
        with open(arg) as f:
            lines = f.readlines()  
        # exclude lines that are not transition matrix     
        lines = lines[-4:]
        minmat = np.empty((len(lines), len(lines)))
        # convert from a text matrix to a numpy array
        for i, line in enumerate(lines): 
            minmat[i, :] = np.squeeze(np.array(pd.read_csv(StringIO(line[:-1]),header=None,  sep=" "))).astype(float)

    # get sequence to compare to
    elif opt in ["-i", "--input"]:        
        for record in SeqIO.parse(arg, "fasta"):
                se = record.seq
                input = list(str(se))   

# get states into an ordered dict
states_dict = OrderedDict()
for i, val in enumerate(states):   
    states_dict[val] = i

# neutral values for addition/multiplication
log_odds = 0
p_plus = 1
p_min = 1
# combine matrices
allmats = np.array([plusmat, minmat])
for i, val in enumerate(input):
    # at the beginning, insert probability from beginning to respective letter
    if i == 0: 
        c = states_dict[val]
        r = states_dict["*"]
    else:
        c = states_dict[val]
        r = states_dict[input[i-1]]
        
    p_p = allmats[0, r,c]
    p_m = allmats[1, r,c]
    p_plus = p_plus * p_p
    p_min = p_min * p_m
    log_odds = log_odds + np.log(p_p/p_m)
 
c = states_dict["+"]
r = states_dict[val]
p_p = allmats[0, r,c]
p_m = allmats[1, r,c]
p_plus = p_plus * p_p
p_min = p_min * p_m
log_odds = log_odds + np.log(p_p/p_m)

print("The probability of the minus model is: ", p_min)
print("The probability of the plus model is: ", p_plus)
print("The log-odds ratio is: ", log_odds)


