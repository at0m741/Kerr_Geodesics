import numpy as np

def compute_riemann_tensor(r, theta, G=1, c=1, m=1):
    riemann = {}
    riemann["R_t_r_t_r"] = 2 * G * m / (c**2 * r**3 - 2 * G * m * r**2)
    riemann["R_t_r_r_t"] = -riemann["R_t_r_t_r"]
    riemann["R_t_theta_t_theta"] = -G * m / (c**2 * r)
    riemann["R_t_theta_theta_t"] = -riemann["R_t_theta_t_theta"]
    riemann["R_t_phi_t_phi"] = -G * m * np.sin(theta)**2 / (c**2 * r)
    riemann["R_t_phi_phi_t"] = -riemann["R_t_phi_t_phi"]
    riemann["R_r_t_t_r"] = (2 * (G * c**2 * m * r - 2 * G**2 * m**2)) / (c**4 * r**4)
    riemann["R_r_t_r_t"] = -riemann["R_r_t_t_r"]
    riemann["R_r_theta_r_theta"] = -G * m / (c**2 * r)
    riemann["R_r_theta_theta_r"] = -riemann["R_r_theta_r_theta"]
    riemann["R_r_phi_r_phi"] = -G * m * np.sin(theta)**2 / (c**2 * r)
    riemann["R_r_phi_phi_r"] = -riemann["R_r_phi_r_phi"]
    riemann["R_theta_t_t_theta"] = -(G * c**2 * m * r - 2 * G**2 * m**2) / (c**4 * r**4)
    riemann["R_theta_t_theta_t"] = -riemann["R_theta_t_t_theta"]
    riemann["R_theta_r_r_theta"] = G * m / (c**2 * r**3 - 2 * G * m * r**2)
    riemann["R_theta_r_theta_r"] = -riemann["R_theta_r_r_theta"]
    riemann["R_theta_phi_theta_phi"] = 2 * G * m * np.sin(theta)**2 / (c**2 * r)
    riemann["R_theta_phi_phi_theta"] = -riemann["R_theta_phi_theta_phi"]
    riemann["R_phi_t_t_phi"] = -(G * c**2 * m * r - 2 * G**2 * m**2) / (c**4 * r**4)
    riemann["R_phi_t_phi_t"] = -riemann["R_phi_t_t_phi"]
    riemann["R_phi_r_r_phi"] = G * m / (c**2 * r**3 - 2 * G * m * r**2)
    riemann["R_phi_r_phi_r"] = -riemann["R_phi_r_r_phi"]
    riemann["R_phi_theta_theta_phi"] = -2 * G * m / (c**2 * r)
    riemann["R_phi_theta_phi_theta"] = -riemann["R_phi_theta_theta_phi"]
    return riemann

r = 20
theta = np.pi / 3.8
riemann_tensor = compute_riemann_tensor(r, theta)

for key, value in riemann_tensor.items():
    print(f"{key}: {value:.6f}")
