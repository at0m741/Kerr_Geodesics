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
