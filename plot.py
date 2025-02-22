import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("simulation_data.csv")
plt.figure(figsize=(8, 6))
plt.subplot(2, 1, 1)
plt.plot(data["time"], data["K_trace"], 'b-')
plt.xlabel("Time")
plt.ylabel("Trace K")
plt.subplot(2, 1, 2)
plt.plot(data["time"], data["Hamiltonian_constraint"], 'r-')
plt.xlabel("Time")
plt.ylabel("Hamiltonian Constraint")
plt.tight_layout()
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_csv("grid_data.csv")

r_values = np.unique(df["r"])
theta_values = np.unique(df["theta"])
R, Theta = np.meshgrid(r_values, theta_values)

components = ["gamma11", "gamma22", "gamma33"]

# Créer une figure avec plusieurs sous-plots
fig, axes = plt.subplots(1, len(components), figsize=(18, 6))

for i, comp in enumerate(components):
    Z_sim = df.pivot(index="theta", columns="r", values=comp).values

    # Formule analytique
    if comp == "gamma33":
        Z_analytic = (R**2) * (np.sin(Theta)**2)
        print(Z_analytic)
    elif comp == "gamma11":
        Z_analytic = np.ones_like(R)  
        print(Z_analytic)
        print(R)
    elif comp == "gamma22":
        Z_analytic = R**2  
        print(Z_analytic)
        print(R)

    # Affichage des données
    ax = axes[i]
    c = ax.pcolormesh(R, Theta, Z_sim, shading='auto', cmap="inferno", alpha=0.6)
    ax.contour(R, Theta, Z_analytic, colors="cyan", linewidths=1)
    
    ax.set_xlabel("r")
    ax.set_ylabel("theta")
    ax.set_title(f"Comparaison de {comp}")

    fig.colorbar(c, ax=ax, label=comp)

# Configuration générale
plt.suptitle("Comparaison des composantes gamma_ij : Simulation vs Théorie")
plt.tight_layout()
plt.show()
