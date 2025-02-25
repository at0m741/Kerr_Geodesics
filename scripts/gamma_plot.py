import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("gamma_slice_full.csv")

NX = int(np.sqrt(len(data)))
NZ = NX

x = data["x"].values.reshape((NX, NZ))
z = data["z"].values.reshape((NX, NZ))

gamma_labels = [
    "gamma_00", "gamma_01", "gamma_02",
    "gamma_10", "gamma_11", "gamma_12",
    "gamma_20", "gamma_21", "gamma_22"
]
gamma_values = [data[label].values.reshape((NX, NZ)) for label in gamma_labels]

fig, axes = plt.subplots(3, 3, figsize=(12, 12))
fig.suptitle("Composantes γ_ij dans une coupe 2D", fontsize=16)

vmin = min(np.min(g) for g in gamma_values)
vmax = max(np.max(g) for g in gamma_values)

for idx, ax in enumerate(axes.flat):
    im = ax.contourf(x, z, gamma_values[idx], levels=100, vmin=vmin, vmax=vmax)
    ax.set_title(gamma_labels[idx])
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    ax.grid(True, linestyle="--", linewidth=0.5, color="white")

fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.8, label="Valeurs γ_ij")

plt.tight_layout()
plt.show()
