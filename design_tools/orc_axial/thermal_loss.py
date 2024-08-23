import math

# Given dummy values
diameter = 0.22  # Diameter of the cylinder in meters
length = 0.17  # Length of the cylinder in meters
T_s = 450.0  # Surface temperature of the cylinder in Kelvin
T_infinity = 310.0  # Ambient air temperature in Kelvin

# Constants
g = 9.81  # Acceleration due to gravity in m/s^2
nu = 15.89e-6  # Kinematic viscosity of air at 300 K in m^2/s
k = 0.0262  # Thermal conductivity of air at 300 K in W/m·K
Pr = 0.71  # Prandtl number for air at 300 K
beta = 1 / T_infinity  # Coefficient of thermal expansion for air (1/K)

# Characteristic length (diameter of the cylinder)
L_c = diameter

# Surface area of the cylinder
A_s = math.pi * diameter * length

# Grashof number (Gr)
Gr = (g * beta * (T_s - T_infinity) * L_c**3) / nu**2

# Rayleigh number (Ra)
Ra = Gr * Pr

# Nusselt number (Nu) using correlation for horizontal cylinder
Nu = 0.36 * (Ra**(1/4))

# Convective heat transfer coefficient (h)
h = (Nu * k) / L_c

# Heat loss due to natural convection (Q_conv)
Q_conv = h * A_s * (T_s - T_infinity)

# Print the results
print(f"Surface area of the cylinder (A_s): {A_s:.4f} m^2")
print(f"Grashof number (Gr): {Gr:.2e}")
print(f"Rayleigh number (Ra): {Ra:.2e}")
print(f"Nusselt number (Nu): {Nu:.2f}")
print(f"Convective heat transfer coefficient (h): {h:.4f} W/m^2·K")
print(f"Heat loss due to natural convection (Q_conv): {Q_conv:.2f} W")
