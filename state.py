import numpy as np

# Créer des grilles pour rho et temp
rho = np.logspace(-8, 8, 100)  # densité en g/cm^3
temp = np.logspace(0, 8, 100)  # température en K

# Créer des valeurs factices pour la pression et l'énergie interne
pressure = np.outer(rho, temp)  # pression en dyne/cm^2 (exemple simplifié)
internal_energy = np.outer(rho, temp**2)  # énergie interne en erg/g (exemple simplifié)

# Sauvegarder dans un fichier npz
np.savez('eos_table.npz', rho=rho, temp=temp, pressure=pressure, internal_energy=internal_energy)

import numpy as np

eos_filename = 'eos_table.npz'
eos_data = np.load(eos_filename)
print(eos_data.files)
