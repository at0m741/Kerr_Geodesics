import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Charger le fichier CSV
file_path = "gamma_slice_full.csv"
df = pd.read_csv(file_path)

# Extraire les coordonnées x, z
x = df["x"].values
z = df["z"].values

# Liste des composantes de γ_ij
gamma_labels = [
    "gamma_00", "gamma_01", "gamma_02",
    "gamma_10", "gamma_11", "gamma_12",
    "gamma_20", "gamma_21", "gamma_22"
]

# Conversion des données en grille
NX = len(np.unique(x))  
NZ = len(np.unique(z))  

# Création de la figure
fig, axes = plt.subplots(3, 3, figsize=(15, 12))  
fig.subplots_adjust(hspace=0.3, wspace=0.3)

for idx, col in enumerate(gamma_labels):
    data = df[col].values.reshape(NX, NZ)
    ax = axes[idx // 3, idx % 3]
    im = ax.imshow(data, extent=[x.min(), x.max(), z.min(), z.max()], origin="lower", aspect="auto", cmap="coolwarm")
    ax.set_title(col)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    fig.colorbar(im, ax=ax, label=col)

plt.suptitle("Composantes de γ_{ij} (Metric Tensor Spatial Components)", fontsize=16)
plt.show()
