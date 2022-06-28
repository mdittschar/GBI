import numpy as np
import sys
import getopt
import Bio
from Bio import SeqIO

argv = sys.argv[1:]
# get a tuple of the arguments to use
opts, _ = getopt.getopt(argv, "a:b:i:", ['matrix1', 'matrix2', 'input'])


for opt, arg in opts:
    if opt in ["-a", "--matrix1"]:
        plusmat = np.load(arg)
    if opt in ["-b", "--matrix2"]:
        minmat = np.load(arg)
    # get new transition matrix name
    elif opt in ["-i", "--input"]:
        
        for record in SeqIO.parse(arg, "fasta"):
                se = record.seq
                input = list(str(se))       


allmats = np.array([plusmat, minmat])
for i, val in enumerate(input):
    if i == 0: 
        if val == "G":
            p_plus = allmats[0, 2, 0]
            p_min = allmats[1, 2, 0]
        elif val == "C":
            p_plus = allmats[0, 2, 1]
            p_min = allmats[1, 2, 1]
        log_odds = np.log(p_plus/p_min)
    elif i < len(input) - 1:
        if val == "G":
            if input[i - 1] == "G":
                p_p = allmats[0, 0, 0]
                p_m = allmats[1, 0, 0]
                p_plus = p_plus * p_p
                p_min = p_min * p_m
            elif input[i - 1] == "C":
                p_p = allmats[0, 1,0]
                p_m = allmats[1, 1, 0]
                p_plus = p_plus * p_p
                p_min = p_min * p_m
            log_odds = log_odds + np.log(p_p/p_m)
        elif val == "C":
            if input[i - 1] == "G":
                p_p = allmats[0, 0, 1]
                p_m = allmats[1, 0, 1]
                p_plus = p_plus * p_p
                p_min = p_min * p_m
            elif input[i - 1] == "C":
                p_p = allmats[0, 1, 1]
                p_m = allmats[1, 1, 1]
                p_plus = p_plus * p_p
                p_min = p_min * p_m
            log_odds = log_odds + np.log(p_p/p_m)
        
    else: 
        if val == "G":
            p_plus = p_plus * allmats[0, 0, 3]
            p_min = p_min * allmats[1, 0, 3]
        elif val == "C":
            p_plus = p_plus * allmats[0, 1, 3]
            p_min = p_min * allmats[1, 1, 3]
        log_odds = log_odds + np.log(allmats[0, 1, 3]/allmats[1, 1, 3])

print("The probability of the minus model is: ", p_min)
print("The probability of the plus model is: ", p_plus)
print("The log-odds ratio is: ", log_odds)

