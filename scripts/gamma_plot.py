import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

file_path = "gamma_slice_full.csv"
df = pd.read_csv(file_path)

x = df["x"].values
z = df["z"].values

gamma_labels = [
    "gamma_00", "gamma_01", "gamma_02",
    "gamma_10", "gamma_11", "gamma_12",
    "gamma_20", "gamma_21", "gamma_22"
]

vals = [df[col].values for col in gamma_labels]
all_vals = np.concatenate(vals)
p_low, p_high = np.percentile(all_vals, [1, 99])
print("Plage dynamique (1-99 percentiles):", p_low, p_high)

NX = len(np.unique(x))
NZ = len(np.unique(z))

fig, axes = plt.subplots(3, 3, figsize=(15, 12))
fig.subplots_adjust(hspace=0.3, wspace=0.3)

for idx, col in enumerate(gamma_labels):
    data = df[col].values.reshape(NX, NZ)
    ax = axes[idx // 3, idx % 3]
    
    center = 1.0 if col in ["gamma_00", "gamma_11", "gamma_22"] else 0.0
    
    norm = mcolors.TwoSlopeNorm(vmin=p_low, vcenter=center, vmax=p_high)
    
    im = ax.imshow(
        data,
        extent=[x.min(), x.max(), z.min(), z.max()],
        origin="lower",
        aspect="auto",
        cmap="viridis",      
        norm=norm,
        interpolation="bilinear" 
    )
    
    ax.set_title(col, fontsize=12)
    ax.set_xlabel("x", fontsize=10)
    ax.set_ylabel("z", fontsize=10)
    fig.colorbar(im, ax=ax, label=col)

plt.suptitle("Composantes de γ_{ij} (Metric Tensor Spatial Components)", fontsize=16)
plt.tight_layout()
plt.show()
