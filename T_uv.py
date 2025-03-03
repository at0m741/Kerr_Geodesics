
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

file_path = "T_energy_momentum.csv"
df = pd.read_csv(file_path)

x = df["x"].values
z = df["z"].values

T_labels = [
    "T_00", "T_01", "T_02",
    "T_10", "T_11", "T_12",
    "T_20", "T_21", "T_22"
]

NX = len(np.unique(x))
NZ = len(np.unique(z))

vals = []
for col in T_labels:
    vals.append(df[col].values)
all_vals = np.concatenate(vals)
vmin, vmax = all_vals.min() / 2, all_vals.max() / 2

fig, axes = plt.subplots(3, 3, figsize=(15, 12))
fig.subplots_adjust(hspace=0.05, wspace=0.15)

for idx, col in enumerate(T_labels):
    data = df[col].values.reshape(NX, NZ)
    ax = axes[idx // 3, idx % 3]
    im = ax.imshow(
        data,
        extent=[x.min(), x.max(), z.min(), z.max()],
        origin="lower",
        aspect="auto",
        cmap="inferno",
        vmin=vmin,    
        vmax=vmax
    )
    ax.set_title(col)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    fig.colorbar(im, ax=ax, label=col)

plt.suptitle("Tenseur énergie-impulsion \( T_{ab} \)", fontsize=16)
plt.tight_layout()
plt.show()
