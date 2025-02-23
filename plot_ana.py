
#!/usr/bin/env python3
import sympy as sp

def main():
    r, theta, phi, M, a = sp.symbols('r theta phi M a', real=True)
    
    Sigma = r**2 + a**2 * sp.cos(theta)**2
    Delta = r**2 - 2*M*r + a**2

    gamma_rr       = Sigma/Delta
    gamma_thetatheta = Sigma
    gamma_phiphi   = (r**2 + a**2 + (2*M*a**2*r*sp.sin(theta)**2)/Sigma) * sp.sin(theta)**2

    gamma = sp.Matrix([
        [gamma_rr,           0,               0],
        [0,          gamma_thetatheta,         0],
        [0,                   0,         gamma_phiphi]
    ])

    print("3-Métrique spatiale gamma_ij :")
    sp.pretty_print(gamma)
    print("\n")

    gamma_inv = gamma.inv()
    print("Inverse de la 3-métrique gamma^ij :")
    sp.pretty_print(gamma_inv)
    print("\n")

    # Liste des coordonnées spatiales
    coords = [r, theta, phi]

    Gamma = {}
    for i in range(3):
        for j in range(3):
            for k in range(3):
                Gamma[i, j, k] = 0
                for l in range(3):
                    term = sp.diff(gamma[l, k], coords[j]) + sp.diff(gamma[l, j], coords[k]) - sp.diff(gamma[j, k], coords[l])
                    Gamma[i, j, k] += gamma_inv[i, l] * term
                Gamma[i, j, k] = sp.simplify(0.5 * Gamma[i, j, k])

    print("Symboles de Christoffel Gamma^i_{jk} de la 3-métrique :\n")
    for i in range(3):
        for j in range(3):
            for k in range(3):
                print(f"Gamma^{i}_{j}{k} =")
                sp.pretty_print(Gamma[i, j, k])
                print("\n")

if __name__ == '__main__':
    main()
