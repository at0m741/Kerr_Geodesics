import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv("gauge_slice.csv")

x_unique = np.sort(df['x'].unique())
z_unique = np.sort(df['z'].unique())
NX = len(x_unique)
NZ = len(z_unique)

alpha = df['alpha'].values.reshape((NX, NZ))
beta0 = df['beta0'].values.reshape((NX, NZ))
beta1 = df['beta1'].values.reshape((NX, NZ))
beta2 = df['beta2'].values.reshape((NX, NZ))
print("alpha: ", alpha)
d_alpha_dt = df['d_alpha_dt'].values.reshape((NX, NZ))
d_beta0_dt = df['d_beta0_dt'].values.reshape((NX, NZ))
d_beta1_dt = df['d_beta1_dt'].values.reshape((NX, NZ))
d_beta2_dt = df['d_beta2_dt'].values.reshape((NX, NZ))
print("d_alpha_dt: ", d_alpha_dt)
fig, axs = plt.subplots(2, 4, figsize=(14, 7))

fields = [(alpha, "Alpha"), (beta0, "Beta_0"), (beta1, "Beta_1"), (beta2, "Beta_2"),
          (d_alpha_dt, "d(Alpha)/dt"), (d_beta0_dt, "d(Beta_0)/dt"), (d_beta1_dt, "d(Beta_1)/dt"), (d_beta2_dt, "d(Beta_2)/dt")]

for ax, (data, title) in zip(axs.flat, fields):
    im = ax.imshow(data, extent=[z_unique.min(), z_unique.max(), x_unique.min(), x_unique.max()],
                   origin='lower', aspect='auto', cmap='coolwarm')
    ax.set_title(title)
    fig.colorbar(im, ax=ax)

plt.tight_layout()
plt.show()
