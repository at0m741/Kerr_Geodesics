
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

filename = "hamiltonian_slice.csv"
df = pd.read_csv(filename)

z_fixed = 0
df_slice = df[df["z"] == z_fixed]

if df_slice.empty:
    print(f"Aucune donnée pour Z = {z_fixed}")
    exit()

x_unique = np.sort(df_slice["x"].unique())
y_unique = np.sort(df_slice["y"].unique())

X, Y = np.meshgrid(x_unique, y_unique)
H = df_slice.pivot(index="y", columns="x", values="Hamiltonian").values

H_min, H_max = np.nanmin(H), np.nanmax(H)
if H_max - H_min == 0:  
    H_norm = np.zeros_like(H)
else:
    H_norm = (H - H_min) / (H_max - H_min)

H_norm = np.nan_to_num(H_norm, nan=0.0)
plt.figure(figsize=(8, 6))
plt.imshow(H, extent=[X.min(), X.max(), Y.min(), Y.max()], origin='lower', cmap="plasma", alpha=0.1 + 1.9 * H_norm)

plt.colorbar(label="Hamiltonian Constraint")

plt.xlabel("X")
plt.ylabel("Y")
plt.title(f"Hamiltonian Constraint Heatmap (Z = {z_fixed}) avec Transparence")

plt.show()
