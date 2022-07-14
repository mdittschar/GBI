#from gff import gffParser
import sys
import getopt
from dittschar_auckenthaler_mygff import myGffParser
import numpy as np

# functions for accuracy, sensitivity and specificity
def acc(tp, tn, fp, fn):
    return (tp + tn)/(tp+tn+fp+fn)

def sensitivity(tp, fn):
    return tp/(tp+fn)

def specificity(tn, fp):
    return tn/(tn + fp)

# get the input arguments
input_file = sys.argv[1:]
opts, _ = getopt.getopt(input_file, "a:b:r:", ['prokka', 'genemark', 'reference'])
for opt, arg in opts:
    if opt in ["-a", "--prokka"]:
        with open(arg) as f:
            lines = f.readlines()
    elif opt in ["-b", "--genemark"]:
        with open(arg) as f:
            lines2 = f.readlines()
    elif opt in ["-r", "--reference"]:
        with open(arg) as f:
            reference = f.readlines()

# parse Prokka
out = myGffParser(lines)
# parse Genemark
out2 = myGffParser(lines2)
# parse reference
ref = myGffParser(reference)

# initialise prediction arrays with all zeros
pred_arr1_neg = np.zeros((6264404))
pred_arr1_pos = np.zeros((6264404))
for i in out.df.index:
    # make sure type is cds
    if out.df.loc[i, "type"] == "CDS":
        # consider both strands
        if out.df.loc[i, "strand"] == "-":
            # mark cds from index start to finish
            pred_arr1_neg[out.df.loc[i, "start"]:out.df.loc[i, "end"]] = 1
        elif out.df.loc[i, "strand"] == "+":
            pred_arr1_pos[out.df.loc[i, "start"]:out.df.loc[i, "end"]] = 1

# initialise prediction arrays
pred_arr2_neg = np.zeros((6264404))
pred_arr2_pos = np.zeros((6264404))
for i in out2.df.index:
    if out2.df.loc[i, "type"] == "CDS":
        # consider both strands
        if out2.df.loc[i, "strand"] == "-":
            # mark cds from index start to finish
            pred_arr2_neg[out2.df.loc[i, "start"]:out2.df.loc[i, "end"]] = 1
        elif out2.df.loc[i, "strand"] == "+":
            pred_arr2_pos[out2.df.loc[i, "start"]:out2.df.loc[i, "end"]] = 1

# initialise prediction arrays
ref_arr_neg = np.zeros((6264404))
ref_arr_pos = np.zeros((6264404))
for i in ref.df.index:
    # make sure type is cds
    if ref.df.loc[i, "type"] == "CDS":
        # consider both strands
        if ref.df.loc[i, "strand"] == "-":
            ref_arr_neg[ref.df.loc[i, "start"]:ref.df.loc[i, "end"]] = 1
        elif ref.df.loc[i, "strand"] == "+":
            ref_arr_pos[ref.df.loc[i, "start"]:ref.df.loc[i, "end"]] = 1



tp_1 = np.sum(np.logical_and(pred_arr1_neg == 1, ref_arr_neg == 1)) + np.sum(np.logical_and(pred_arr1_pos == 1, ref_arr_pos == 1))#/np.sum(ref_arr)
print("True positive rate of PROKKA: ", tp_1)

tn_1 = np.sum(np.logical_and(pred_arr1_neg == 0, ref_arr_neg == 0)) + np.sum(np.logical_and(pred_arr1_pos == 0, ref_arr_pos == 0))#/(len(ref_arr) - np.sum(ref_arr))
print("True negative rate of PROKKA: ", tn_1)

fp_1 = np.sum(np.logical_and(pred_arr1_neg == 1, ref_arr_neg == 0)) + np.sum(np.logical_and(pred_arr1_pos == 1, ref_arr_pos == 0))#/(len(ref_arr) - np.sum(ref_arr))
print("False positive rate of PROKKA: ", fp_1)#

fn_1 = np.sum(np.logical_and(pred_arr1_neg == 0, ref_arr_neg == 1)) + np.sum(np.logical_and(pred_arr1_pos == 0, ref_arr_pos == 1))#/np.sum(ref_arr)
print("False negative rate of PROKKA: ", fn_1)

acc_1 = acc(tp_1, tn_1, fp_1, fn_1)
print("The accuracy of PROKKA is :", acc_1)

sens_1 = sensitivity(tp_1, fn_1)
print("The sensitivity of PROKKA is: ", sens_1)

spec_1 = specificity(tn_1, fp_1)
print("The specificity of PROKKA is: ", spec_1)
print(f"Total length is: {tp_1+tn_1+fp_1+fn_1}")
print()

# values for Genemark output
tp_2 = np.sum(np.logical_and(pred_arr2_neg == 1, ref_arr_neg == 1)) + np.sum(np.logical_and(pred_arr2_pos == 1, ref_arr_pos == 1))#/np.sum(ref_arr)
print("True positive rate of Genemark: ", tp_2)

tn_2 = np.sum(np.logical_and(pred_arr2_neg == 0, ref_arr_neg == 0)) + np.sum(np.logical_and(pred_arr2_pos == 0, ref_arr_pos == 0))#/(len(ref_arr) - np.sum(ref_arr))
print("True negative rate of Genemark: ", tn_2)

fp_2 = np.sum(np.logical_and(pred_arr2_neg == 1, ref_arr_neg == 0)) + np.sum(np.logical_and(pred_arr2_pos == 1, ref_arr_pos == 0))#/(len(ref_arr) - np.sum(ref_arr))
print("False positive rate of Genemark: ", fp_2)#
fn_2 = np.sum(np.logical_and(pred_arr2_neg == 0, ref_arr_neg == 1))+np.sum( np.logical_and(pred_arr2_pos == 0, ref_arr_pos == 1))#/np.sum(ref_arr)
print("False negative rate of Genemark: ", fn_2)

acc_2 = acc(tp_2, tn_2, fp_2, fn_2)
print("The accuracy of Genemark is :", acc_2)

sens_2 = sensitivity(tp_2, fn_2)
print("The sensitivity of Genemark is: ", sens_2)

spec_2 = specificity(tn_2, fp_2)
print("The specificity of Genemark is: ", spec_2)


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