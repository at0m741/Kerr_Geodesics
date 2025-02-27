
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Charger le fichier CSV
file_path = "gauge_slice.csv"
df = pd.read_csv(file_path)

# Extraire les coordonnées x, z
x = df["x"].values
z = df["z"].values

# Liste des champs de jauge et dérivées temporelles
gauge_labels = [
    "alpha", "beta0", "beta1", "beta2",
    "d_alpha_dt", "d_beta0_dt", "d_beta1_dt", "d_beta2_dt"
]

# Conversion des données en grille
NX = len(np.unique(x))  
NZ = len(np.unique(z))  

# Création de la figure
fig, axes = plt.subplots(2, 4, figsize=(18, 8))  
fig.subplots_adjust(hspace=0.3, wspace=0.3)

for idx, col in enumerate(gauge_labels):
    data = df[col].values.reshape(NX, NZ)
    ax = axes[idx // 4, idx % 4]
    im = ax.imshow(data, extent=[x.min(), x.max(), z.min(), z.max()], origin="lower", aspect="auto", cmap="coolwarm")
    ax.set_title(col)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    fig.colorbar(im, ax=ax, label=col)

plt.suptitle("Champs de jauge et dérivées temporelles", fontsize=16)
plt.show()
