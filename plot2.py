
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_christoffel_component(df, k, a, b):
    """
    Visualise la composante Gamma3[k][a][b] en fonction de (r, theta).
    k, a, b sont des entiers entre 0 et 2.
    """
    col_name = f"Gamma3_{k}_{a}_{b}"

    # Vérifier que la colonne existe
    if col_name not in df.columns:
        print(f"La colonne {col_name} n'existe pas dans le DataFrame.")
        return

    pivot_df = df.pivot(index="theta", columns="r", values=col_name)

    R_values = pivot_df.columns.values
    theta_values = pivot_df.index.values
    Z = pivot_df.values

    # Construire la grille pour contourf
    R_grid, Theta_grid = np.meshgrid(R_values, theta_values)

    plt.figure(figsize=(8, 6))
    # Affichage par contour
    cont = plt.contourf(R_grid, Theta_grid, Z, levels=50, cmap="viridis")
    plt.colorbar(cont, label=col_name)
    plt.xlabel("r")
    plt.ylabel("theta")
    plt.title(f"Visualisation de {col_name} (r, theta)")
    plt.show()

def main():
    # Charger le fichier CSV
    df = pd.read_csv("christoffel_data.csv")

    print(df.head())

    plot_christoffel_component(df, k=0, a=1, b=1)
    plot_christoffel_component(df, k=0, a=2, b=2)
    plot_christoffel_component(df, k=1, a=0, b=0)
    # Tu peux visualiser d'autres composantes :
    # plot_christoffel_component(df, k=1, a=1, b=1)
    # plot_christoffel_component(df, k=2, a=2, b=0)
    # etc.

if __name__ == "__main__":
    main()
