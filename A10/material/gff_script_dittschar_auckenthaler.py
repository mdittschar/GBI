from gff import gffParser
import sys
import getopt
from mygff import myGffParser
import numpy as np


out = myGffParser()
print(out.df["start"][1000])

pred_arr = np.zeros((len(out.df.index)))
for i in out.df.index:
    pred_arr[out.df["start"][i]:out.df["end"][i]] = 1
print("Pred_arr: ", pred_arr)

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