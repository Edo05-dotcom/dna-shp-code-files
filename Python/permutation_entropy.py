import numpy as np
import antropy as ant

# File stores time series of transcrition events
FILE = "output/history.dat" 

ORDER = 5 # length of permutation
DELAY = 1
NORMALIZE = True # divide by log(M!)
WINDOW = 20

data = np.loadtxt(FILE)
gene_ids = data[:, 1].astype(int)

genes = np.unique(gene_ids)

pe_list = []
weights = []

for g in genes:
    # Build transcription amount time series for gene g
    nwin = len(gene_ids) // WINDOW

    activity = np.array([
        np.sum(gene_ids[i*WINDOW:(i+1)*WINDOW] == g)
        for i in range(nwin)
    ], dtype=float)

    pe = ant.perm_entropy(activity,
                          order=ORDER,
                          delay=DELAY,
                          normalize=NORMALIZE)

    pe_list.append(pe)
    weights.append(np.sum(activity))  # total activity of gene

# Weighted average
pe_weighted = np.average(pe_list, weights=weights)

print("Weighted permutation entropy =", pe_weighted)