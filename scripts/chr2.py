import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Charger le fichier CSV
file_path = "christoffel_slice.csv"
df = pd.read_csv(file_path)

# Extraire les coordonnées x, z
x = df["x"].values
z = df["z"].values

# Liste des symboles de Christoffel à conserver
christoffel_labels = [
    "Gamma000", "Gamma001", "Gamma002",
    "Gamma010", "Gamma011", 
    "Gamma020", "Gamma021", "Gamma022",
    "Gamma100", "Gamma101",
    "Gamma110", "Gamma111", "Gamma112",
    "Gamma120", "Gamma121",
    "Gamma200", "Gamma202",
    "Gamma211", "Gamma212",
    "Gamma220", "Gamma221", "Gamma222"
]

# Conversion des données en grille
NX = len(np.unique(x))  
NZ = len(np.unique(z))  

# Création de la figure
fig, axes = plt.subplots(3, 7, figsize=(20, 12))  
fig.subplots_adjust(hspace=0.3, wspace=0.3)

for idx, col in enumerate(christoffel_labels):
    data = df[col].values.reshape(NX, NZ)
    ax = axes[idx // 7, idx % 7]
    im = ax.imshow(data, extent=[x.min(), x.max(), z.min(), z.max()], origin="lower", aspect="auto", cmap="coolwarm")
    ax.set_title(col)
    ax.set_xlabel("x")
    ax.set_ylabel("z")
    fig.colorbar(im, ax=ax, label=col)

plt.suptitle("Symboles de Christoffel (Sans Gamma210, Gamma201, Gamma012, Gamma102, Gamma120)", fontsize=16)
plt.show()
