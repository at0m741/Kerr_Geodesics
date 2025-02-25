
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Charger les données
data = pd.read_csv("alpha_slice.csv")

# Extraire NX et NZ dynamiquement
NX = len(data["x"].unique())
NZ = len(data["z"].unique())

# Reconstruction des grilles
x = data["x"].values.reshape((NX, NZ))
z = data["z"].values.reshape((NX, NZ))
alpha = data["alpha"].values.reshape((NX, NZ))

# Création de la figure
fig, ax = plt.subplots(figsize=(8, 6))
fig.suptitle("Lapse \(\\alpha\) dans une coupe 2D", fontsize=16)

# Normalisation des couleurs pour voir les variations
vmin, vmax = np.percentile(alpha, [2, 98])

c = ax.pcolormesh(x, z, alpha, shading='auto', cmap='plasma', vmin=vmin, vmax=vmax)
ax.set_xlabel("x")
ax.set_ylabel("z")
fig.colorbar(c, ax=ax, label="Valeur de \(\\alpha\)")

plt.tight_layout()
plt.show()
