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

vals = []
for col in gamma_labels:
    vals.append(df[col].values)
all_vals = np.concatenate(vals)
vmin, vmax = all_vals.min(), all_vals.max()

mean_val = all_vals.mean()
delta = 0.0001  # +/- 1%
vmin = mean_val - delta
vmax = mean_val + delta
NX = len(np.unique(x))
NZ = len(np.unique(z))

fig, axes = plt.subplots(3, 3, figsize=(15, 12))
fig.subplots_adjust(hspace=0.05, wspace=0.15)

for idx, col in enumerate(gamma_labels):
    data = df[col].values.reshape(NX, NZ)
    ax = axes[idx // 3, idx % 3]
    im = ax.imshow(
        data,
        extent=[x.min(), x.max(), z.min(), z.max()],
        origin="lower",
        aspect="auto",
        cmap="coolwarm",
        vmin=vmin,    
        vmax=vmax
    )
    ax.set_title(col)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    fig.colorbar(im, ax=ax, label=col)

plt.suptitle("Composantes de γ_{ij} (Metric Tensor Spatial Components)", fontsize=16)
plt.tight_layout()
plt.show()
