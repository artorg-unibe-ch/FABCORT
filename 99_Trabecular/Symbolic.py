import sympy as sp

# Identity tensor
I = sp.eye(3)

# Green-Lagrange strain tensor in 3D space
F11, F12, F13 = sp.symbols('F_{11} F_{12} F_{13}', positive=True)
F21, F22, F23 = sp.symbols('F_{21} F_{22} F_{23}', positive=True)
F31, F32, F33 = sp.symbols('F_{31} F_{32} F_{33}', positive=True)

F = sp.Matrix([[F11, F12, F13],
               [F21, F22, F23],
               [F31, F32, F33]])

E = 1/2 * (F + F.T) - I

# Cauchy stress tensor in 3D space
T11, T12, T13 = sp.symbols('T_{11} T_{12} T_{13}', positive=True)
T21, T22, T23 = sp.symbols('T_{21} T_{22} T_{23}', positive=True)
T31, T32, T33 = sp.symbols('T_{31} T_{32} T_{33}', positive=True)

T = sp.Matrix([[T11, T12, T13],
               [T21, T22, T23],
               [T31, T32, T33]])

# Rewrite tensors in 6D space
E6 = sp.Matrix([E[0,0],
                E[1,1],
                E[2,2],
                2*E[1,2],
                2*E[0,2],
                2*E[0,1]])

T6 = sp.Matrix([T[0,0],
                T[1,1],
                T[2,2],
                T[1,2],
                T[0,2],
                T[0,1]])

# Elasticity tensor in 3D space
C = sp.zeros(6,6)
for i in range(6):
    for j in range(6):
        C[i,j] = T6[i] / E6[j]

# Convert components into 6D space
C6 = sp.zeros(6,6)
for i in range(3):
    for j in range(3):
        C6[i,j] = C[i,j]

for i in range(3):
    for j in range(3):
        C6[i+3,j+3] = 2 * C[i+3,j+3]

for i in range(3):
    for j in range(3,6):
        C6[i,j] = 2**0.5 * C[i,j]

for i in range(3,6):
    for j in range(3):
        C6[i,j] = 2**0.5 * C[i,j]

# Orthotropic symmetry
for i in range(6):
    for j in range(6):
        if (i >= 3 or j >= 3) and i!=j:
            C6[i,j] = 0

# Eigenvalue and eigenvectors
Basis = C6.eigenvects()