import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import math

# Charger les données depuis le CSV
data = pd.read_csv("christoffel_slice.csv")

all_cols = data.columns.tolist()
columns_to_plot = all_cols[2:]  

vals = [data[col].values for col in columns_to_plot]
all_vals = np.concatenate(vals)
p_low, p_high = np.percentile(all_vals, [1, 99])
print("Plage dynamique (1er et 99ème percentiles):", p_low, p_high)

nplots = len(columns_to_plot)

for i in range(0, nplots, 2):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    for j in range(2):
        idx = i + j
        if idx >= nplots:
            axes[j].axis("off")
            continue
            
        col = columns_to_plot[idx]
        pivoted = data.pivot(index="z", columns="x", values=col)
        
        ax = axes[j]
        norm = mcolors.TwoSlopeNorm(vmin=p_low, vcenter=0, vmax=p_high)
        
        im = ax.imshow(
            pivoted.values,
            extent=[pivoted.columns.min(), pivoted.columns.max(), pivoted.index.min(), pivoted.index.max()],
            origin="lower",
            cmap='inferno',
            norm=norm,
            interpolation='bilinear'
        )
        
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
    
    plt.tight_layout()
    plt.show()
