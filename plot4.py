import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplcursors

# Charger le fichier CSV
csv_file = "output/extrinsic_curvature.csv"
df = pd.read_csv(csv_file)

# Vérification des valeurs
print(df.head())

# Grouper par (i, j)
grouped = df.groupby(['i', 'j'])
tensor_dict = {}

for (i, j), group in grouped:
    tensor_str = ""
    for a in range(3):
        lignes = []
        for b in range(3):
            val_array = group[(group['a'] == a) & (group['b'] == b)]['Kij'].values
            if len(val_array) > 0:
                lignes.append(f"K[{a}][{b}] = {val_array[0]:.2e}")
            else:
                lignes.append(f"K[{a}][{b}] = NaN")
        tensor_str += "\n".join(lignes) + "\n\n"
    tensor_dict[(i, j)] = tensor_str

# Création de la figure
fig, ax = plt.subplots(figsize=(10, 8))
ax.set_title("Cliquez ou survolez un point pour afficher le tenseur K_{ij}")
ax.set_xlabel("j (coordonnée en θ)")
ax.set_ylabel("i (coordonnée en r)")
ax.grid(True)

# Extraction des valeurs uniques pour les axes
i_vals = sorted(df['i'].unique())
j_vals = sorted(df['j'].unique())
ax.set_xticks(j_vals)
ax.set_yticks(i_vals)

# Tracé des points
points = np.array([[j, i] for (i, j) in tensor_dict.keys()])
sc = ax.scatter(points[:, 0], points[:, 1], s=100, c='blue')

# Ajout du curseur interactif
cursor = mplcursors.cursor(sc, hover=True)

@cursor.connect("add")
def on_add(sel):
    x, y = sel.target
    j = int(round(x))
    i = int(round(y))
    tensor_str = tensor_dict.get((i, j), "Aucun tenseur trouvé")
    sel.annotation.set(text=tensor_str, fontsize=8)

plt.tight_layout()
plt.show()

# Calcul des normes L2 de K_{ij}
df_kij = df[['i', 'j', 'Kij']].groupby(['i', 'j']).sum().reset_index()
norm_dict = {(row['i'], row['j']): np.sqrt(row['Kij']**2) for _, row in df_kij.iterrows()}

# Création du DataFrame pour le heatmap
norm_df = pd.DataFrame({
    'i': [i for (i, j) in norm_dict.keys()],
    'j': [j for (i, j) in norm_dict.keys()],
    'norm': [v for v in norm_dict.values()]
})
heatmap_norm = norm_df.pivot(index='i', columns='j', values='norm')

# Affichage du heatmap
fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
im = ax.imshow(heatmap_norm, origin='lower', aspect='auto', cmap='viridis',
               extent=[min(j_vals)-0.5, max(j_vals)+0.5, min(i_vals)-0.5, max(i_vals)+0.5],
               interpolation='bicubic')

ax.set_title("Heatmap : Norme L2 du tenseur K_{ij} par cellule")
ax.set_xlabel("j (coordonnée en θ)")
ax.set_ylabel("i (coordonnée en r)")
ax.set_xticks(j_vals)
ax.set_yticks(i_vals)
fig.colorbar(im, ax=ax, label="Norme L2")

plt.tight_layout()
plt.show()
