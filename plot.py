#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def animate_gamma111(frame, scat, ax):
    # Construire le nom du fichier pour le frame donné
    filename = f"output/grid_t{frame:04d}.csv"
    if not os.path.exists(filename):
        print(f"File {filename} not found.")
        return scat

    try:
        df = pd.read_csv(filename)
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return scat

    # On suppose que le CSV contient les colonnes 'r', 'theta' et 'Gamma3_1_1_1'
    # Ces colonnes doivent être générées lors de l'enregistrement de la grille
    r = df['r'].values
    theta = df['theta'].values
    gamma111 = df['Gamma11'].values

    # Met à jour les positions et la couleur des points dans le scatter
    offsets = np.column_stack((theta, r))
    scat.set_offsets(offsets)
    scat.set_array(gamma111)
    ax.set_title(f"Gamma3[1][1][1] at frame {frame}")
    return scat

def main():
    initial_filename = "output/grid_t0000.csv"
    try:
        df0 = pd.read_csv(initial_filename)
    except Exception as e:
        print(f"Error reading {initial_filename}: {e}")
        return

    # Extraction des coordonnées et de la valeur initiale de Gamma3_1_1_1
    theta0 = df0['theta'].values
    r0 = df0['r'].values
    gamma111_0 = df0['Gamma3_1_1_1'].values

    # Création de la figure et du scatter plot
    fig, ax = plt.subplots(figsize=(8, 6))
    scat = ax.scatter(theta0, r0, c=gamma111_0, cmap='viridis', s=20)
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel('r')
    ax.set_title("Gamma3[1][1][1] at frame 0")
    cbar = fig.colorbar(scat, ax=ax, label=r'$\Gamma^1_{11}$ value')
    
    # Nombre de frames à animer
    n_frames = 30
    ani = animation.FuncAnimation(
        fig,
        animate_gamma111,
        frames=n_frames,
        fargs=(scat, ax),
        interval=500,  # intervalle en ms entre frames
        repeat=True
    )

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
