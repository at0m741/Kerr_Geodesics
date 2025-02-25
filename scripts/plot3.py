import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("K_slice.csv")

x_unique = np.sort(df['x'].unique())
z_unique = np.sort(df['z'].unique())
NX = len(x_unique)
NZ = len(z_unique)

fig, axs = plt.subplots(3, 3, figsize=(12, 10))

components = [['K00', 'K01', 'K02'],
              ['K10', 'K11', 'K12'],
              ['K20', 'K21', 'K22']]

for i in range(3):
    for j in range(3):
        comp = components[i][j]
        data = df[comp].values.reshape((NX, NZ))
        im = axs[i, j].imshow(data, extent=[z_unique.min(), z_unique.max(), x_unique.min(), x_unique.max()],
                                origin='lower', aspect='auto', cmap='viridis')
        axs[i, j].set_title(comp)
        fig.colorbar(im, ax=axs[i, j])
        
plt.tight_layout()
plt.show()
