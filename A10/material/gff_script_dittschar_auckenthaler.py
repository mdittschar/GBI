from gff import gffParser
import sys
import getopt
from mygff import myGffParser
import numpy as np

def acc(tp, tn, fp, fn):
    return (tp + tn)/(tp+tn+fp+fn)

def sensitivity(tp, fn):
    return tp/(tp+fn)

def specificity(tn, fp):
    return tn/(tn + fp)


input_file = sys.argv[1:]
opts, _ = getopt.getopt(input_file, "a:b:r:", ['input1', 'input2', 'reference'])
for opt, arg in opts:
    if opt in ["-a", "--input1"]:
        with open(arg) as f:
            lines = f.readlines()
    elif opt in ["-b", "--input2"]:
        with open(arg) as f:
            lines2 = f.readlines()
    elif opt in ["-r", "--reference"]:
        with open(arg) as f:
            reference = f.readlines()

out = myGffParser(lines)
out2 = myGffParser(lines2)
ref = myGffParser(reference)
print(ref.df)
#print(out.df["start"][1000])

pred_arr1 = np.zeros((6264404))
for i in out.df.index:
    pred_arr1[out.df["start"][i]:out.df["end"][i]] = 1

pred_arr2 = np.zeros((6264404))
for i in out.df.index:
    pred_arr2[out2.df["start"][i]:out2.df["end"][i]] = 1

ref_arr = np.zeros((6264404))
for i in out.df.index:
    if ref.df["type"][i] == "CDS":
        ref_arr[ref.df["start"][i]:ref.df["end"][i]] = 1
print("Pred_arr: ", ref_arr)


tp_1 = np.mean(np.logical_and(pred_arr1 == 1, ref_arr == 1))
print("True positive rate of prokka: ", tp_1)

tn_1 = np.mean(np.logical_and(pred_arr1 == 0, ref_arr == 0))
print("True negative rate of prokka: ", tn_1)

fp_1 = np.mean(np.logical_and(pred_arr1 == 1, ref_arr == 0))
print("False positive rate of prokka: ", fp_1)#

fn_1 = np.mean(np.logical_and(pred_arr1 == 0, ref_arr == 1))
print("False negative rate of prokka: ", fn_1)

acc_1 = acc(tp_1, fn_1, fp_1, fn_1)
print("The accuracy of prokka is :", acc_1)

sens_1 = sensitivity(tp_1, fn_1)
print("The sensitivity of PROKKA is: ", sens_1)

spec_1 = specificity(tn_1, fp_1)
print("The specificity of PROKKA is: ", spec_1)


# input_file = sys.argv[1:]
# opts, _ = getopt.getopt(input_file, "i:", ['input'])
# print(opts)
# for opt, arg in opts:
#     if opt in ["-i", "--input"]:
#         with open(arg) as f:
#             lines = f.readlines()
                    
# out = gffParser(lines)
# print(out.data["NC_002516.2"][0]["strand"])#["frame"])#[out.data["NC_002516.2"][0]["feature"] == "CDS"])
# print(out.data["NC_002516.2"][10]["feature"])#["feature"] == "CDS")
#print(out.getGenes('NC_002516.2'))

# ## get genes in the chromosome 1
#print(out.data["MPBCFPNP_00001"])

# ## get mRNA corresponding to a gene
# out.getmRNA(chrom, gene)

# ## get coding regions in the mRNA
# out.getCDS(chrom, mRNA)

# ## get introns/exons
# out.getInrons(chrom, mRNA)
# out.getExonss(chrom, mRNA) 