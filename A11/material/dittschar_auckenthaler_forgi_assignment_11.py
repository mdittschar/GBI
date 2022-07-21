import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi
cg = forgi.load_rna("auckenthaler_dittschar_sequence1.fx", allow_many=False)
fvm.plot_rna(cg, text_kwargs={"fontweight":"black"}, lighten=0.7,
             backbone_kwargs={"linewidth":3})
plt.show()

cg = forgi.load_rna("auckenthaler_dittschar_sequence2.fx", allow_many=False)
fvm.plot_rna(cg, text_kwargs={"fontweight":"black"}, lighten=0.7,
             backbone_kwargs={"linewidth":3})
plt.show()