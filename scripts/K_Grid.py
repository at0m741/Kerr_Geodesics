import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("K_full.csv")
x_vals = np.unique(df["x"])
y_vals = np.unique(df["y"])
z_vals = np.unique(df["z"])

NX, NY, NZ = len(x_vals), len(y_vals), len(z_vals)

K00 = df["K00"].values.reshape(NX, NY, NZ)
K11 = df["K11"].values.reshape(NX, NY, NZ)
K22 = df["K22"].values.reshape(NX, NY, NZ)
K_trace = K00 + K11 + K22

fig, axes = plt.subplots(1, 5, figsize=(20,5))

yslices = np.linspace(0, NY-1, 5, dtype=int)

for i, y_idx in enumerate(yslices):
    im = axes[i].imshow(K_trace[:, y_idx, :], 
                        origin='lower', 
                        extent=[x_vals[0], x_vals[-1], z_vals[0], z_vals[-1]],
                        cmap='viridis', aspect='auto')
    axes[i].set_title(f"y={y_vals[y_idx]:.2f}")
    axes[i].set_xlabel("x")
    axes[i].set_ylabel("z")
fig.colorbar(im, ax=axes.ravel().tolist(), label="K trace")
plt.show()
