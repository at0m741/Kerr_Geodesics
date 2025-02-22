import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("grid_data.csv")

M = 1.0  
a = 0.9  

def Sigma(r, theta):
    return r**2 + a**2 * np.cos(theta)**2

def Delta(r):
    return r**2 - 2 * M * r + a**2

r_values = np.unique(df["r"])
theta_values = np.unique(df["theta"])
R, Theta = np.meshgrid(r_values, theta_values)

components = [
    "gamma11", "gamma22", "gamma33",
    "Gamma_1_0_0", "Gamma_1_1_1", "Gamma_1_2_2", "Gamma_2_2_1" 
]

n_rows = 2
n_cols = (len(components) + 1) // 2 

fig, axes = plt.subplots(n_rows, n_cols, figsize=(24, 12))
axes = axes.reshape(n_rows, n_cols)

comparison_results = []

for i, comp in enumerate(components):
    Z_sim = df.pivot(index="theta", columns="r", values=comp).values

    if comp == "gamma33":
        Z_analytic = (R**2 + a**2 + (2 * M * a**2 * R * np.sin(Theta)**2) / Sigma(R, Theta)) * np.sin(Theta)**2
    elif comp == "gamma11":
        Z_analytic = Sigma(R, Theta) / Delta(R)
    elif comp == "gamma22":
        Z_analytic = Sigma(R, Theta)
    elif comp == "Gamma_1_0_0":
        Z_analytic = (M * (2 * R - a**2 * np.sin(Theta)**2)) / (Sigma(R, Theta) * Delta(R))
    elif comp == "Gamma_1_1_1":
        Z_analytic = (M * (R**2 - a**2 * np.cos(Theta)**2)) / (Sigma(R, Theta) * Delta(R))
    elif comp == "Gamma_1_2_2":
        Z_analytic = - Delta(R) / Sigma(R, Theta)
    elif comp == "Gamma_2_2_1":
        Z_analytic = (a**2 * np.sin(Theta) * np.cos(Theta) * (R**2 + a**2 + 2 * M * a**2 * R / Sigma(R, Theta))) / Sigma(R, Theta)

    row, col = divmod(i, n_cols)
    ax = axes[row, col]
    c = ax.pcolormesh(R, Theta, Z_sim, shading='auto', cmap="inferno", alpha=0.6)
    ax.contour(R, Theta, Z_analytic, colors="cyan", linewidths=1)

    max_diff = np.max(np.abs(Z_sim - Z_analytic))
    min_sim = np.min(Z_sim)
    max_sim = np.max(Z_sim)
    min_analytic = np.min(Z_analytic)
    max_analytic = np.max(Z_analytic)

    comparison_results.append({
        "Composante": comp,
        "Max Différence": max_diff,
        "Min Simulée": min_sim,
        "Max Simulée": max_sim,
        "Min Analytique": min_analytic,
        "Max Analytique": max_analytic
    })

    ax.set_xlabel("r")
    ax.set_ylabel("theta")
    ax.set_title(f"Comparaison de {comp} (Kerr)")

    fig.colorbar(c, ax=ax, label=comp)

plt.suptitle("Comparaison des composantes de la métrique et des symboles de Christoffel : Simulation vs Théorie (Métrique de Kerr)")
plt.tight_layout()
plt.show()

comparison_df = pd.DataFrame(comparison_results)
print(comparison_df)
