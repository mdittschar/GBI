from gff import gffParser
import sys
    
input_file = sys.argv[1]
out = gffParser(input_file)

# ## get genes in the chromosome 1
#print(out.data["MPBCFPNP_00001"])

# ## get mRNA corresponding to a gene
# out.getmRNA(chrom, gene)

# ## get coding regions in the mRNA
# out.getCDS(chrom, mRNA)

# ## get introns/exons
# out.getInrons(chrom, mRNA)
# out.getExonss(chrom, mRNA) 