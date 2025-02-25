import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("christoffel_slice_full.csv")

NX = int(np.sqrt(len(data)))
NZ = NX

x = data['x'].values.reshape((NX, NZ))
z = data['z'].values.reshape((NX, NZ))

gamma_cols = []
for kk in range(3):
    for aa in range(3):
        for bb in range(3):
            col_name = f"Gamma{kk}_{aa}{bb}"
            gamma_cols.append(col_name)

for col_name in gamma_cols:
    arr = data[col_name].values.reshape((NX, NZ))

    plt.figure(figsize=(8,6))
    plt.contourf(x, z, arr, levels=60, cmap="inferno")
    plt.colorbar(label=col_name)
    plt.xlabel("x")
    plt.ylabel("z")
    plt.title(f"Christoffel slice : {col_name}")
    plt.show()
