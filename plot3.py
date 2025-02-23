import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import mplcursors

csv_file = "output/christoffel.csv"
df = pd.read_csv(csv_file, header=None, skiprows=1, names=['i', 'j', 'k', 'l', 'm', 'Gamma3', 'Ricci', 'Kij'])

grouped = df.groupby(['i', 'j'])
tensor_dict = {}
for (i, j), group in grouped:
    tensor_str = ""
    for k in range(3):
        lignes = []
        for l in range(3):
            valeurs = []
            for m in range(3):
                val_array = group[(group['k'] == k) & (group['l'] == l) & (group['m'] == m)]['Gamma3'].values 
                val_array2 = group[(group['k'] == k) & (group['l'] == l) & (group['m'] == m)]['Ricci'].values
                val_array3 = group[(group['k'] == k) & (group['l'] == l) & (group['m'] == m)]['Kij'].values
                if len(val_array) > 0:
                    valeurs.append(f"{val_array[0]:.2e}")
                else:
                    valeurs.append("NaN")
                if len(val_array2) > 0:
                    valeurs.append(f"{val_array2[0]:.2e}")
                else:
                    valeurs.append("NaN")
            lignes.append("[" + ", ".join(valeurs) + "]")
        tensor_str += f"k={k}:\n" + "\n".join(lignes) + "\n\n"
    tensor_dict[(i, j)] = tensor_str

fig, ax = plt.subplots(figsize=(10, 8))
ax.set_title("Cliquez ou survolez un point pour afficher le tenseur Γ³")
ax.set_xlabel("j (coordonnée en θ)")
ax.set_ylabel("i (coordonnée en r)")
ax.grid(True)

i_vals = sorted(df['i'].unique())
j_vals = sorted(df['j'].unique())
ax.set_xticks(j_vals)
ax.set_yticks(i_vals)

points = np.array([[j, i] for (i, j) in tensor_dict.keys()])
sc = ax.scatter(points[:, 0], points[:, 1], s=100, c='blue')

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

df_gamma = df[['i', 'j', 'Gamma3']].groupby(['i', 'j']).sum().reset_index()
df_ricci = df[['i', 'j', 'Ricci']].groupby(['i', 'j']).sum().reset_index()
df_kij = df[['i', 'j', 'Kij']].groupby(['i', 'j']).sum().reset_index()
norm_dict = {(row['i'], row['j']): np.sqrt(row['Gamma3']**2) for _, row in df_gamma.iterrows()}
norm_dict2 = {(row['i'], row['j']): np.sqrt(row['Ricci']**2) for _, row in df_ricci.iterrows()}
norm_dict3 = {(row['i'], row['j']): np.sqrt(row['Kij']**2) for _, row in df_kij.iterrows()}
norm_df = pd.DataFrame({'i': [i for (i, j) in norm_dict.keys()], 'j': [j for (i, j) in norm_dict.keys()], 'norm': [v for v in norm_dict.values()]})
heatmap_norm = norm_df.pivot(index='i', columns='j', values='norm')

norm_df2 = pd.DataFrame({'i': [i for (i, j) in norm_dict2.keys()], 'j': [j for (i, j) in norm_dict2.keys()], 'norm': [v for v in norm_dict2.values()]})
heatmap_norm2 = norm_df2.pivot(index='i', columns='j', values='norm')

norm_df3 = pd.DataFrame({'i': [i for (i, j) in norm_dict3.keys()], 'j': [j for (i, j) in norm_dict3.keys()], 'norm': [v for v in norm_dict3.values()]})
heatmap_norm3 = norm_df3.pivot(index='i', columns='j', values='norm')


fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
im = ax.imshow(heatmap_norm, origin='lower', aspect='auto', cmap='viridis',
               extent=[min(j_vals)-0.5, max(j_vals)+0.5, min(i_vals)-0.5, max(i_vals)+0.5],
               interpolation='bicubic')

ax.set_title("Heatmap : Norme L2 du tenseur Γ³ par cellule")
ax.set_xlabel("j (coordonnée en θ)")
ax.set_ylabel("i (coordonnée en r)")
ax.set_xticks(j_vals)
ax.set_yticks(i_vals)
fig.colorbar(im, ax=ax, label="Norme L2")
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
im2 = ax.imshow(heatmap_norm2, origin='lower', aspect='auto', cmap='plasma',
                extent=[min(j_vals)-0.5, max(j_vals)+0.5, min(i_vals)-0.5, max(i_vals)+0.5],
                interpolation='bicubic')

ax.set_title("Heatmap : Norme L2 du tenseur Ricci par cellule")
ax.set_xlabel("j (coordonnée en θ)")
ax.set_ylabel("i (coordonnée en r)")
ax.set_xticks(j_vals)
ax.set_yticks(i_vals)
fig.colorbar(im2, ax=ax, label="Norme L2")
plt.tight_layout()
plt.show()


fig, ax = plt.subplots(figsize=(8, 6), dpi=150)
im3 = ax.imshow(heatmap_norm3, origin='lower', aspect='auto', cmap='plasma',
                extent=[min(j_vals)-0.5, max(j_vals)+0.5, min(i_vals)-0.5, max(i_vals)+0.5],
                interpolation='bicubic')
ax.set_title("Heatmap : Norme L2 du tenseur Kij par cellule")
ax.set_xlabel("j (coordonnée en θ)")
ax.set_ylabel("i (coordonnée en r)")
ax.set_xticks(j_vals)
ax.set_yticks(i_vals)
fig.colorbar(im3, ax=ax, label="Norme L2")

plt.tight_layout()
plt.show()

