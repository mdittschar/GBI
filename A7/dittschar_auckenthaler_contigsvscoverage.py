import matplotlib.pyplot as plt
import numpy as np

c = np.arange(0,20, .1)
thetas = [.2,.3,.4,.5,.6]

for theta in thetas:
    N = np.exp(-c*(1-theta)) *c
    plt.plot(c,np.exp(-c*(1-theta))*N, label=f"theta = {theta}")
plt.title("Expected Number of Contigs by Coverage and Minimal Overlap")
plt.xlabel("Coverage")
plt.ylabel("Expected number of contigs")
plt.legend()
plt.show()