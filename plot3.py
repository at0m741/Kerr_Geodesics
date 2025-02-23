import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Charger les données
filename = "output/dK_evolution.csv"
df = pd.read_csv(filename)

if "t" not in df.columns or "i" not in df.columns or "j" not in df.columns or "r" not in df.columns or "theta" not in df.columns or "dK_11" not in df.columns:
    raise ValueError("Le fichier CSV doit contenir les colonnes : t, i, j, r, theta, dK_22")

Nr = df["i"].max() + 1
Ntheta = df["j"].max() + 1

timesteps = sorted(df["t"].unique())

frames = []
for t in timesteps:
    frame = np.full((Nr, Ntheta), np.nan)  
    subset = df[df["t"] == t]
    for _, row in subset.iterrows():
        frame[int(row["i"]), int(row["j"])] = row["dK_01"]
        frame[int(row["i"]), int(row["j"])] = row["dK_11"]
    frames.append(frame)

# Initialiser la figure
fig, ax = plt.subplots(figsize=(8, 6))
cmap = plt.get_cmap("inferno")  # Palette de couleurs
im = ax.imshow(frames[0], cmap=cmap, origin="lower", extent=[df["theta"].min(), df["theta"].max(), df["r"].min(), df["r"].max()])
ax.set_xlabel("Theta")
ax.set_ylabel("r")
ax.set_title("Évolution de dK_11")

# Ajouter une barre de couleur
cbar = plt.colorbar(im, ax=ax)
cbar.set_label("dK_11")

# Fonction pour mettre à jour l'animation
def update(frame_index):
    im.set_array(frames[frame_index])
    ax.set_title(f"Évolution de dK_11 - t = {timesteps[frame_index]}")
    return [im]

# Créer l'animation
ani = animation.FuncAnimation(fig, update, frames=len(frames), interval=100, blit=False)

# Afficher l'animation
plt.show()
