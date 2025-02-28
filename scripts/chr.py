import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

# Lecture du CSV
data = pd.read_csv("christoffel_slice.csv")

all_cols = data.columns.tolist() 
columns_to_plot = all_cols[2:] 

nplots = len(columns_to_plot) 
ncols = 3
nrows = math.ceil(nplots / ncols)

fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 5*nrows), squeeze=False)

for idx, col in enumerate(columns_to_plot):
    pivoted = data.pivot(index="z", columns="x", values=col)
    xvals = pivoted.columns.values
    zvals = pivoted.index.values
    X, Z = np.meshgrid(xvals, zvals)

    ax = axes[idx // ncols, idx % ncols]
    im = ax.pcolormesh(X, Z, pivoted.values, shading='auto', cmap='RdBu')
    ax.set_title(col)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    fig.colorbar(im, ax=ax, label=col)

for idx_ax in range(nplots, nrows*ncols):
    ax = axes[idx_ax // ncols, idx_ax % ncols]
    ax.set_visible(False)

plt.tight_layout()
plt.show()
