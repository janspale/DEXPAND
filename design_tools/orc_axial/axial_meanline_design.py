"""
This is a Python design tool for the design of an axial impulse single stage small scale ORC turbine.
Accompanies the PhD thesis of J.Spale, 2024.
Jan.Spale@cvut.cz
"""
from __future__ import print_function

# Importing the libraries
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import timeit
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
from tabulate import tabulate
from units import convert

# Set the REFPROP path
RP_path = "c:\\Program Files (x86)\\REFPROP"

# cycle boundary conditions
p_inlet = convert(550,"kPa","Pa") # kPa
T_inlet = convert(180,"C","K") # C
p_outlet = convert(55,"kPa","Pa") # kPa
m_dot = 0.25  # kg/s
fluid = "MM" # hexamethyldisiloxane
z = [1.0] # [-] molar fraction of the fluid

# Design parameters
n_design = 18000  # rpm
D_mid = 0.135  # m
alpha_stator = 13  # deg
h_over_D_mid = 0.04  # [-]

#initialize REFPROP
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX']) # Instantiate the REFPROP function library
RP.SETPATHdll(os.environ['RPPREFIX']) # Set the path to the folder containing the REFPROP shared library
RP.SETFLUIDSdll(fluid)  # Set the fluids
baseSI = RP.GETENUMdll(0, "MASSBASESI").iEnum # Get the base SI units
iMass = 0; iFlag = 0 # Set the mass and flag to 0

def check_err(r): # Function to check the error code from REFPROP
    if r.ierr > 0:
        raise ValueError(r.herr)

#initialize state variables
h,s,T,p,rho,a,Z,Ma = np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5)

# turbine inlet ([0])
p[0] = p_inlet # Pa
T[0] = T_inlet # K
r = RP.REFPROPdll(fluid,"PT","H;S;D;W;Z",baseSI,iMass,iFlag,p[0],T[0],z); check_err(r)
h[0],s[0],rho[0],a[0] = r.Output[0:4] # J/kg, J/kgK, kg/m3, m/s
Ma[0] = 0 # [-], Mach number at the inlet

# stator isentropic expansion ([1])
p[1] = p_outlet # Pa
s[1] = s[0] # J/kgK, isentropic entropy
r = RP.REFPROPdll(fluid,"PS","H;T;D;W;Z",baseSI,iMass,iFlag,p[1],s[1],z); check_err(r) # isentropic expansion
h[1],T[1],rho[1],a[1] = r.Output[0:4] # J/kg, K, kg/m3, m/s

#velocities
dh_stage = h[0]-h[1] # J/kg, isentropic enthalpy drop
c_is = np.sqrt(2*dh_stage) # m/s, isentropic velocity
U = n_design*np.pi*D_mid/60 # m/s, mean blade speed midspan
U_over_c_is = U/c_is # [-], ratio of mean blade speed to isentropic velocity

Ma[1] = c_is/a[1]  # [-], Mach number at the stator outlet (isentropic)
phi_stator = np.sqrt(1-(0.0029*Ma[1]**3-0.0502*Ma[1]**2+0.2241*Ma[1]-0.0877)) # [-], velocity loss coefficient at the nozzles

""" obsolete Mach number iteration
# sonic point (nozzle throat) ([2])
#get isentropic exponent for (0) based on T,rho
r = RP.REFPROPdll(fluid,"TD","CP/CV",baseSI,iMass,iFlag,T[0],rho[0],z); check_err(r)
gamma_star = r.Output[0] # [-], isentropic exponent
T_star_is = T[0]*(2/(gamma_star+1)) # K, isentropic temperature at the nozzle throat
T_star = T_star_is # K, temperature at the nozzle throat, iterated, initial guess is isentropic

#get speed of sound, density and enthalpy and pressure at the nozzle throat
r = RP.REFPROPdll(fluid,"TD","D;H;W;P",baseSI,iMass,iFlag,T_star,rho[0],z); check_err(r)
rho[2],h[2],a[2],p[2] = r.Output[0:4] # kg/m3, J/kg, m/s, Pa
"""
""" obsolete T_star iteration
#iterate to find the actual position of the nozzle throat
it = 0
loop_cond = False
Cond1 = False
Cond2 = False

# Placeholders for variables that need initial values

step = 0.01  # Adjust as necessary
max_it = 10000  # Adjust as necessary

while not loop_cond:
    it += 1
    #rho as a function of T_star and s[0]
    r = RP.REFPROPdll(fluid,"TS","D",baseSI,iMass,iFlag,T_star,s[0],z); check_err(r)
    rho[2] = r.Output[0] # kg/m3

    #h,a,p as a function of T_star and rho[2]
    r = RP.REFPROPdll(fluid,"TD","H;W;P",baseSI,iMass,iFlag,T_star,rho[2],z); check_err(r)
    h[2],a[2],p[2] = r.Output[0:3] # kg/m3, J/kg, m/s, Pa

    # finding the value of T_star that yields such a value of h[2] that a[2] is equal to c_is
    X = (h[0] - h[2])
    Y = (a[2] ** 2) / 2

    comparison_value = (X - Y) / Y # [-], 
    if comparison_value > 0: # T_star is too low
        T_star += step
        Cond1 = True
    else: # T_star is too high
        T_star -= step
        Cond2 = True

    if it >= max_it or (Cond1 and Cond2): # if the iteration limit is reached or the conditions are met
        loop_cond = True

print("Number of iterations: ", it)
print("Initial guess: ", T_star_is)
print("T_star: ", T_star)
"""

#calculate isentropic expansion quasi-steady 1D, 1000 steps in pressure
p_is = np.linspace(p[0],p[1],1000)
T_is,rho_is,h_is,a_is,c_is,Ma_is,PR_is,v_is = np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is))

for i in range(len(p_is)):
    r = RP.REFPROPdll(fluid,"PS","T;D;H;W",baseSI,iMass,iFlag,p_is[i],s[0],z); check_err(r)
    T_is[i],rho_is[i],h_is[i],a_is[i] = r.Output[0:4] # K, kg/m3, J/kg, m/s
    c_is[i] = np.sqrt(2*(h[0]-h_is[i])) # m/s
    Ma_is[i] = c_is[i]/a_is[i] # [-]
    PR_is[i] = p_is[i]/p_outlet # [-]
    v_is[i] = 1/rho_is[i] # m3/kg

#locate properties at the nozzle throat where Ma = 1
Ma_is_1 = np.where(Ma_is > 1)[0][0]
T_is_1 = T_is[Ma_is_1]
p_is_1 = p_is[Ma_is_1]
rho_is_1 = rho_is[Ma_is_1]
h_is_1 = h_is[Ma_is_1]
a_is_1 = a_is[Ma_is_1]
c_is_1 = c_is[Ma_is_1]

#calculate the actual expansion quasi-steady 1D, 1000 steps in pressure
#nozzle out enthalpy
h_nozzle_out = h[0] - ((phi_stator**2)*(h[0]-h_is[-1]))# J/kg

#get nozzle out T, rho, a,s, Ma,c
r = RP.REFPROPdll(fluid,"HP","T;D;W;S",baseSI,iMass,iFlag,h_nozzle_out,p_outlet,z); check_err(r)
T_nozzle_out,rho_nozzle_out,a_nozzle_out,s_nozzle_out = r.Output[0:4] # K, kg/m3, m/s
c_nozzle_out = phi_stator*c_is[-1] # m/s
Ma_nozzle_out = c_nozzle_out/a_nozzle_out # [-]

#TODO: draw expansion for isentropic case and nozzles with losses case 

""" nozzle with losses vs isentropic expansion 
#calculate the actual expansion quasi-steady 1D using s_nozzle_out and p_act in 1000 steps
p_act = np.linspace(p[0],p_outlet,1000)
s_act = np.zeros(len(p_act)) # J/kgK, actual entropy
#linear production of entropy
T_act,rho_act,h_act,a_act,c_act,Ma_act,PR_act = np.zeros(len(p_act)),np.zeros(len(p_act)),np.zeros(len(p_act)),np.zeros(len(p_act)),np.zeros(len(p_act)),np.zeros(len(p_act)),np.zeros(len(p_act))
for i in range(len(p_act)):
    #if Ma_nozzle_out < 1, s_act = s[0]; if Ma_nozzle_out > 1, s_act linearly increases from s[0] to s_nozzle_out
    r = RP.REFPROPdll(fluid,"PS","T;D;H;W",baseSI,iMass,iFlag,p_act[i],s_nozzle_out,z); check_err(r)
    T_act[i],rho_act[i],h_act[i],a_act[i] = r.Output[0:4] # K, kg/m3, J/kg, m/s
    if h_act[i] > h[0]:
        r=RP.REFPROPdll(fluid,"PS","H",baseSI,iMass,iFlag,p_act[i],s[0],z); check_err(r)
        h_act[i] = r.Output[0] # J/kg
        p_throat = p_act[i]
    else: #linearly interpolate s_act between s[0] and s_nozzle_out p that corresponds to the trigger of the else statement to p_outlet
        s_act[i] = s[0] + (s_nozzle_out-s[0])*(p_act[i]-p_throat)/(p_outlet-p_throat)
        r = RP.REFPROPdll(fluid,"PS","T;D;H;W",baseSI,iMass,iFlag,p_act[i],s_act[i],z); check_err(r)
        T_act[i],rho_act[i],h_act[i],a_act[i] = r.Output[0:4]
    c_act[i] = np.sqrt(2*(h[0]-h_act[i])) # m/s
    Ma_act[i] = c_act[i]/a_act[i] # [-]
    PR_act[i] = p_act[i]/p_outlet # [-]

#locate properties at the nozzle outlet where Ma = 1
Ma_act_1 = np.where(Ma_act > 1)[0][0]
T_act_1 = T_act[Ma_act_1]
p_act_1 = p_act[Ma_act_1]
rho_act_1 = rho_act[Ma_act_1]
h_act_1 = h_act[Ma_act_1]
a_act_1 = a_act[Ma_act_1]
c_act_1 = c_act[Ma_act_1]

#draw isentropic expansion
fig, ax = plt.subplots()
ax.plot(PR_is,Ma_is, label='Isentropic expansion', color='blue')
ax.set(xlabel='Pressure ratio [-]', ylabel='Mach number [-]',
       title='Isentropic expansion')
ax.scatter(p_is_1/p_outlet,1, color='red', label='Throat_is')
ax.plot(PR_act,Ma_act, label='Actual expansion', color='green')
ax.scatter(p_act_1/p_outlet,1, color='orange', label='Throat_act')
ax.grid()
#invert x axis
ax.invert_xaxis()
plt.show()
"""

#calculate nozzle area, throat, outlet and between
A_outlet = m_dot/(rho_nozzle_out*c_nozzle_out) # m2
A_throat = m_dot/(rho_is_1*c_is_1) # m2
A_ratio = A_outlet/A_throat # [-]

#nozzle design
h_nozzle = h_over_D_mid*D_mid # m
b_nozzles_out = A_outlet/h_nozzle # m
b_nozzles_throat = A_throat/h_nozzle # m
no_nozzles = 10 # [-]
b_nozzle_out = b_nozzles_out/no_nozzles # m
b_nozzle_throat = b_nozzles_throat/no_nozzles # m

#==> #design nozzle using open-moc

""" #TODO: how to get Gamma for inlet? Γ=1+(v/a)*(da/dv)s=konst
v1=1/rho[0] # m3/kg
dv = 1e-7 # m3/kg
v2 = v1+dv # m3/kg
a1=a[0] # m/s
a2=RP.REFPROPdll(fluid,"SD","W",baseSI,iMass,iFlag,s[0],1/v2,z).Output[0] # m/s
Gamma = 1 + (v1/a1)*((a2-a1)/dv) # [-]
print("Gamma: ",Gamma)
"""



#compare nozzle design with real gas properties and ideal gas properties

#calculate the velocity triangle nozzle outlet

#calculate the velocity triangle rotor inlet

#calculate the velocity triangle rotor outlet

#calculate rotor loss factor phi_rotor

#draw expansion line

#calculate additional losses

#calculate turbine power

#calculate turbine efficiency

#draw velocity trangles

#run sensitivity analyses

#run optimization using genetic algorithm to find iptimal H/D, alpha,n,D




