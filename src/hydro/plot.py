import numpy as np
import matplotlib.pyplot as plt

# Lire les résultats de la simulation
data = np.loadtxt('output.txt')

# Extraire les variables
rho = data[:, 0]
vx = data[:, 1]
vy = data[:, 2]
vz = data[:, 3]
P = data[:, 4]

#3D plot

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(rho, vx, vy, P, cmap='viridis')
plt.show()