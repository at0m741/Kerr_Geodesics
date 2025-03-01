import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

data = pd.read_csv("christoffel_slice.csv")

all_cols = data.columns.tolist()
columns_to_plot = all_cols[2:]  
vals = []
for col in columns_to_plot:
    vals.append(data[col].values)
all_vals = np.concatenate(vals)
vmin, vmax = all_vals.min(), all_vals.max()

mean_val = all_vals.mean()
delta = 0.1  # +/- 1%
vmin = mean_val - delta
vmax = mean_val + delta
nplots = len(columns_to_plot)
ncols = 2  
nrows = math.ceil(nplots / ncols)

fig, axes = plt.subplots(nrows, ncols, figsize=(7*ncols, 5*nrows))

for idx, col in enumerate(columns_to_plot):
    pivoted = data.pivot(index="z", columns="x", values=col)
    xvals = pivoted.columns.values
    zvals = pivoted.index.values
    X, Z = np.meshgrid(xvals, zvals)

    ax = axes[idx // ncols, idx % ncols]
    im = ax.pcolormesh(X, Z, pivoted.values,
                       shading='auto', cmap='RdBu',
                       vmin=vmin, vmax=vmax)
    M = 1.0 
    r_h = M * (1 + M / 4)**2

    theta = np.linspace(0, 2*np.pi, 300)
    horizon_x = r_h * np.cos(theta)
    horizon_z = r_h * np.sin(theta)

    ax.plot(horizon_x, horizon_z, 'k--', linewidth=2, label="Horizon des événements")
    ax.legend()

    ax.set_title(col)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    fig.colorbar(im, ax=ax, label=col)

for idx_ax in range(nplots, nrows*ncols):
    ax = axes[idx_ax // ncols, idx_ax % ncols]
    ax.set_visible(False)

plt.tight_layout()
plt.show()
