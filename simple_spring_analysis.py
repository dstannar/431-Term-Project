# simple spring model for sage's loading structure
#
# in, lbf

import numpy as np

# Constants
n_legs = 8

height = 10.0             # in
leg_length = 11.15        # in
top_diameter = 3.65       # in

P = 5.0                   # lbf, applied static load / dropped weight
drop_height = 10.0        # in

k_leg = 10.0              # lbf/in, spring coeff of one leg


# angle of leg
cos_theta = height / leg_length


# Stiffness matrix
# One vertical DOF: vertical motion of rigid top disk
K = np.zeros((1, 1))

# Each inclined axial spring contributes k_leg*cos(theta)^2 in vertical stiffness
K[0, 0] = n_legs * k_leg * cos_theta**2


# Static load vector
F = np.zeros(1)
F[0] = P


# Static solution
W = np.linalg.solve(K, F)

static_max_deflection = W[0]
static_Fpeak = K[0, 0] * static_max_deflection


# Dynamic
# Energy balance:
#   0.5*k*delta^2 = P*(drop_height + delta)
#
# P is the dropped weight
# Positive root:
#   delta = (P + sqrt(P^2 + 2*k*P*h))/k

k_total = K[0, 0]

dynamic_max_deflection = (P + np.sqrt(P**2 + 2.0 * k_total * P * drop_height)) / k_total

dynamic_Fpeak = k_total * dynamic_max_deflection


# Print results
print(f"Total vertical stiffness = {k_total:.4f} lbf/in")

print("Static loading:")
print(f"Max deflection = {static_max_deflection:.4f} in")
print(f"Fpeak = {static_Fpeak:.4f} lbf")

print("Dynamic drop loading:")
print(f"Max deflection = {dynamic_max_deflection:.4f} in")
print(f"Fpeak = {dynamic_Fpeak:.4f} lbf")