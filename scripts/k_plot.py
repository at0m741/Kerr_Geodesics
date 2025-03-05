import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

file_path = "K_slice.csv"
df = pd.read_csv(file_path)

x = df["x"].values
z = df["z"].values

K_labels = [
    "K00", "K01", "K02",
    "K10", "K11", "K12",
    "K20", "K21", "K22"
]
print(df.columns)
NX = len(np.unique(x))
NZ = len(np.unique(z))

vals = np.concatenate([df[col].values for col in K_labels])
p_low, p_high = np.percentile(vals, [1, 99])

fig, axes = plt.subplots(3, 3, figsize=(15, 12))
fig.subplots_adjust(hspace=0.3, wspace=0.3)

for idx, col in enumerate(K_labels):
    data = df[col].values.reshape(NX, NZ)
    ax = axes[idx // 3, idx % 3]
    
    center = 0.0
    
    norm = mcolors.TwoSlopeNorm(vmin=p_low, vcenter=center, vmax=p_high)
    
    im = ax.imshow(
        data,
        extent=[x.min(), x.max(), z.min(), z.max()],
        origin="lower",
        aspect="auto",
        cmap="RdBu",
        interpolation="bilinear"
    )
    
    ax.set_title(col, fontsize=12)
    ax.set_xlabel("x", fontsize=10)
    ax.set_ylabel("z", fontsize=10)
    fig.colorbar(im, ax=ax, label=col)

plt.suptitle("Composantes de $K_{ij}$ (Extrinsic Curvature)", fontsize=16)
plt.tight_layout()
plt.show()
