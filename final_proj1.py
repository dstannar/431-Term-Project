# AERO431 Final Structures Project
# FEM analysis for one beam/spring member
# Then combines two stacked springs in series
# and multiple paths in parallel

import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Constants
# ------------------------------------------------------------

E = 3600e6       # Pa, change to your material
nu = 0.35
G = E / (2 * (1 + nu))

#upper spring constant
k_top=1141     #N-m

#lower spring constant
k_bot=1369.2   #N-m

# Geometry of ONE small spring/beam member
L = 0.127       # m, CHANGE to length of one spring
c = 0.0381     # m, width
t = 0.002032     # m, thickness

A = c * t
I = c * t**3 / 12.0
J = (1/3) * c * t**3

# Load
P_total = 5.0 * 4.44822   # 5 lb in N

# Number of identical load paths around structure
n_paths = 8            

# Load per path
P_path = P_total / n_paths

# Since each path has two springs in series,
# analyze ONE spring using the same path load
Px = -P_path
Pz = 0.0
M_end = 0.0

# Mesh
n_elements = 20
x_mesh = np.linspace(0, L, n_elements + 1)
h = L / n_elements

# ------------------------------------------------------------
# 2D frame element stiffness
# local DOFs = [u1, w1, theta1, u2, w2, theta2]
# ------------------------------------------------------------

def frame_element_stiffness(E, A, I, h):
    k = np.zeros((6, 6))

    # Axial stiffness EA
    k[0, 0] = E*A/h
    k[0, 3] = -E*A/h
    k[3, 0] = -E*A/h
    k[3, 3] = E*A/h

    # Bending stiffness EI
    k[1, 1] = 12*E*I/h**3
    k[1, 2] = 6*E*I/h**2
    k[1, 4] = -12*E*I/h**3
    k[1, 5] = 6*E*I/h**2

    k[2, 1] = 6*E*I/h**2
    k[2, 2] = 4*E*I/h
    k[2, 4] = -6*E*I/h**2
    k[2, 5] = 2*E*I/h

    k[4, 1] = -12*E*I/h**3
    k[4, 2] = -6*E*I/h**2
    k[4, 4] = 12*E*I/h**3
    k[4, 5] = -6*E*I/h**2

    k[5, 1] = 6*E*I/h**2
    k[5, 2] = 2*E*I/h
    k[5, 4] = -6*E*I/h**2
    k[5, 5] = 4*E*I/h

    return k

# ------------------------------------------------------------
# Assemble global stiffness matrix for ONE spring
# ------------------------------------------------------------

n_nodes = n_elements + 1
n_dof = 3 * n_nodes

K = np.zeros((n_dof, n_dof))
F = np.zeros(n_dof)

for e in range(n_elements):
    ke = frame_element_stiffness(E, A, I, h)

    dofs = [
        3*e,
        3*e + 1,
        3*e + 2,
        3*(e+1),
        3*(e+1) + 1,
        3*(e+1) + 2
    ]

    for i in range(6):
        for j in range(6):
            K[dofs[i], dofs[j]] += ke[i, j]

# Apply load at free end of one spring
end_node = n_nodes - 1

F[3*end_node] += Px
F[3*end_node + 1] += Pz
F[3*end_node + 2] += M_end

# Fixed left end
fixed_dofs = [0, 1, 2]
free_dofs = np.setdiff1d(np.arange(n_dof), fixed_dofs)

K_reduced = K[np.ix_(free_dofs, free_dofs)]
F_reduced = F[free_dofs]

U = np.zeros(n_dof)
U[free_dofs] = np.linalg.solve(K_reduced, F_reduced)

# ------------------------------------------------------------
# One spring results
# ------------------------------------------------------------

u_tip = U[3*end_node]
w_tip = U[3*end_node + 1]
theta_tip = U[3*end_node + 2]

delta_spring = abs(w_tip)

k_spring = P_path / delta_spring

# ------------------------------------------------------------
# Two stacked springs in series
# ------------------------------------------------------------

delta_path = 2.0 * delta_spring
k_path = k_spring / 2.0

# ------------------------------------------------------------
# Multiple paths in parallel
# ------------------------------------------------------------

k_global = n_paths * k_path
delta_global = P_total / k_global

# ------------------------------------------------------------
# Print results
# ------------------------------------------------------------

print("ONE SPRING FEM RESULTS")
print(f"Load per path = {P_path:.4f} N")
print(f"One spring vertical deflection = {1000*delta_spring:.4f} mm")
print(f"One spring stiffness = {k_spring:.4f} N/m")
print()

print("TWO STACKED SPRINGS IN SERIES")
print(f"Deflection of one full path = {1000*delta_path:.4f} mm")
print(f"Stiffness of one full path = {k_path:.4f} N/m")
print()

print("GLOBAL STRUCTURE")
print(f"Number of paths = {n_paths}")
print(f"Global stiffness = {k_global:.4f} N/m")
print(f"Global deflection under 5 lb = {1000*delta_global:.4f} mm")
print(f"Global deflection under 5 lb = {delta_global/0.0254:.4f} in")

# ------------------------------------------------------------
# Plot one spring deformation
# ------------------------------------------------------------

xs = x_mesh
us = U[0::3]

plt.figure()
plt.plot(xs, 1000*us, 'bo-', lw=2)
plt.xlabel('x (m)')
plt.ylabel('w (mm)')
plt.title('Deflection of One Spring/Beam Member')
plt.grid(True)
plt.show()

#The finite element model shown here represents one spring path as an equivalent beam 
# with stiffness EI. The resulting deflection plot illustrates the overall flexibility 
# of the spring path under load. In the actual structure, deformation would be distributed 
# among the curved segments of the spring, causing the S-shaped geometry to partially 
# straighten while remaining in the elastic range.