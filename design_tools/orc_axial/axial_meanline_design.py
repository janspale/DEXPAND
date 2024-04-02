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
import CoolProp.CoolProp as CP
from tabulate import tabulate
from units import convert
from scipy.interpolate import griddata

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
throat_trigger = False
p_is = np.linspace(p[0],p[1],1000)
T_is,rho_is,h_is,a_is,c_is,Ma_is,PR_is,v_is = np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is))
T_act,rho_act,h_act,a_act,c_act,Ma_act,PR_act,s_act = np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is))
for i in range(len(p_is)):
    r = RP.REFPROPdll(fluid,"PS","T;D;H;W",baseSI,iMass,iFlag,p_is[i],s[0],z); check_err(r)
    T_is[i],rho_is[i],h_is[i],a_is[i] = r.Output[0:4] # K, kg/m3, J/kg, m/s
    c_is[i] = np.sqrt(2*(h[0]-h_is[i])) # m/s
    Ma_is[i] = c_is[i]/a_is[i] # [-]
    PR_is[i] = p_is[i]/p_outlet # [-]
    v_is[i] = 1/rho_is[i] # m3/kg
    #as soon as Ma_is[i] > 1, save the index of the throat as throat_index = i
    if Ma_is[i] > 1 and not throat_trigger:
        throat_index = i
        throat_trigger = True
    if Ma_is[i] < 1:
        h_act[i] = h_is[i] # J/kg
    else:
        h_act[i] = h_is[throat_index] - ((0.875**2)*(h_is[throat_index]-h_is[i])) # J/kg
    r = RP.REFPROPdll(fluid,"PH","T;D;S;W",baseSI,iMass,iFlag,p_is[i],h_act[i],z); check_err(r)
    T_act[i],rho_act[i],s_act[i],a_act[i] = r.Output[0:4] # K, kg/m3, J/kg, m/s
    c_act[i] = np.sqrt(2*(h[0]-h_act[i])) # m/s
    Ma_act[i] = c_act[i]/a_act[i] # [-]
    PR_act[i] = p_is[i]/p_outlet # [-]

#calculate the actual expansion quasi-steady 1D, 1000 steps in pressure
#nozzle out enthalpy
h_nozzle_out = h[0] - ((phi_stator**2)*(h[0]-h_is[-1]))# J/kg

#get nozzle out T, rho, a,s, Ma,c
r = RP.REFPROPdll(fluid,"HP","T;D;W;S",baseSI,iMass,iFlag,h_nozzle_out,p_outlet,z); check_err(r)
T_nozzle_out,rho_nozzle_out,a_nozzle_out,s_nozzle_out = r.Output[0:4] # K, kg/m3, m/s
c_nozzle_out = phi_stator*c_is[-1] # m/s
Ma_nozzle_out = c_nozzle_out/a_nozzle_out # [-]

#TODO: draw expansion for isentropic case and nozzles with losses case 
"""

#draw isentropic expansion
fig, ax = plt.subplots()
ax.plot(PR_is,Ma_is, label='Isentropic expansion', color='black')
ax.set(xlabel='Pressure ratio [-]', ylabel='Mach number [-]', title='Isentropic expansion vs Nozzle with losses')
ax.scatter(p_is[throat_index]/p_outlet,1, color='red', marker='x')
#add text to the throat
ax.text(p_is[throat_index]/p_outlet+0.05, 1.05, 'Throat', ha='right')
ax.plot(PR_act,Ma_act, label='Actual expansion', color='black', linestyle='--')
#add scatter to the nozzle outlet for both and add text with the Ma value, use x markers in red
ax.scatter(p_is[-1]/p_outlet,Ma_act[-1], color='red', marker='x')
ax.text(p_is[-1]/p_outlet, Ma_act[-1], 'Ma = '+str(round(Ma_act[-1],3)), ha='right')
ax.scatter(p_is[-1]/p_outlet,Ma_is[-1], color='red', marker='x')
ax.text(p_is[-1]/p_outlet, Ma_is[-1], 'Ma_is = '+str(round(Ma_is[-1],3)), ha='right')
#add a tick at x axis for the inlet, throat and for the outlet of the nozzle
ax.axvline(p_is[0]/p_outlet, color='black', linestyle='--', linewidth=0.5)
ax.axvline(p_is[throat_index]/p_outlet, color='black', linestyle='--', linewidth=0.5)
ax.axvline(p_is[-1]/p_outlet, color='black', linestyle='--', linewidth=0.5)
#add labels to the ticks, rotate the text by 90°, offset the text by 0.5 in x axis
ax.text(p_is[0]/p_outlet+0.05, 0.5, 'Inlet', ha='right', rotation=90)
ax.text(p_is[throat_index]/p_outlet+0.05, 0.5, 'Throat', ha='right', rotation=90)
ax.text(p_is[-1]/p_outlet+0.05, 0.5, 'Outlet', ha='right', rotation=90)
#write the x axis value next to the tick text, rotated by 90°
ax.text(p_is[0]/p_outlet-0.05, 0.5, "PR = "+str(round(p_is[0]/p_outlet,2)), ha='left', rotation=90)
ax.text(p_is[throat_index]/p_outlet-0.05, 0.5, "PR = "+str(round(p_is[throat_index]/p_outlet,2)), ha='left', rotation=90)
ax.text(p_is[-1]/p_outlet-0.05, 0.5, "PR = "+str(round(p_is[-1]/p_outlet,2)), ha='left', rotation=90)
ax.legend()
ax.grid(True, linestyle='--')
#invert x axis
ax.invert_xaxis()
plt.savefig('isentropic_expansion_nozzle_losses.png', dpi = 600)
#plt.show()
plt.close()

#add ideal gas expansion to the chart for MM?

#calculate nozzle area, throat, outlet and between
A_outlet = m_dot/(rho_nozzle_out*c_nozzle_out) # m2
A_throat = m_dot/(rho_is[throat_index]*c_is[throat_index]) # m2
A_ratio = A_outlet/A_throat # [-]

#nozzle design
h_nozzle = h_over_D_mid*D_mid # m
b_nozzles_out = A_outlet/h_nozzle # m
b_nozzles_throat = A_throat/h_nozzle # m
no_nozzles = 10 # [-]
b_nozzle_out = b_nozzles_out/no_nozzles # m
b_nozzle_throat = b_nozzles_throat/no_nozzles # m
"""

#==> #design nozzle using open-moc

#TODO: how to get Gamma for inlet? Γ=1+(v/a)*(da/dv)s=konst

CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, os.getenv('RPPREFIX'))

"""
# Define the temperature and entropy ranges for the saturated vapor region
T_min, T_max = 325, 560  # Temperature in K, modify as needed
s_min, s_max = 350, 1100 # Entropy in J/(kg*K), modify as needed
steps = 1000  # Number of steps in the grid

# Set up a grid for temperature and entropy
T_gamma = np.linspace(T_min, T_max,steps)  
s_gamma = np.linspace(s_min, s_max,steps)  

# Create a meshgrid for T and s
TT, ss = np.meshgrid(T_gamma, s_gamma)
#strip certain combination of T and s from the meshgrid

# Initialize an array to hold the Gamma values
Gamma = np.zeros_like(TT)
Gamma_line = np.zeros_like(TT)

#get quality to find if the point is in the saturated vapor region
x = np.zeros_like(TT)
#if quality is <= 1, read 0 for Gamma, if quality is > 1, calculate Gamma
for i in range(TT.shape[0]):
    for j in range(TT.shape[1]):
        # Set the state using T and s
        AS = CP.AbstractState("REFPROP", "MM")
        T_check = 0.15*ss[i, j]+410 
        if TT[i, j] < T_check:
            AS.update(CP.SmassT_INPUTS,ss[i, j],TT[i, j])
        else:
            Gamma[i, j] = np.nan
        x[i, j] = AS.Q()
        if x[i, j] <= 1:
            Gamma[i, j] = np.nan
        else:
            Gamma[i, j] = AS.fundamental_derivative_of_gas_dynamics()
        if Gamma[i, j] > 1:
            Gamma[i, j] = np.nan
#read through the meshgrid and find first points where either NaN switches to a number or vice versa and switch that nan to 1
for i in range(Gamma.shape[0]):
    for j in range(Gamma.shape[1]):
        if np.isnan(Gamma[i, j]) and not np.isnan(Gamma[i-1, j]):
            Gamma_line[i, j] = 1
        if not np.isnan(Gamma[i, j]) and np.isnan(Gamma[i-1, j]):
            Gamma_line[i, j] = 1
#calculate saturated vapor line
s_vaps = []
s_liqs = []
Tcrits = []
for Ts in T_gamma:
    if Ts < 518.5:
        Tcrits.append(Ts)
        AS = CP.AbstractState("REFPROP", "MM")
        AS.update(CP.QT_INPUTS,1, Ts)
        s_vapor = AS.smass() # J/(kg*K)
        AS.update(CP.QT_INPUTS,0, Ts)
        s_liquid = AS.smass()
        s_vaps.append(s_vapor)
        s_liqs.append(s_liquid)

#calculate an isobar for p = p_inlet and p = p_outlet
isobar_inlet = []
isobar_outlet = []
for entropy in s_gamma:
    AS = CP.AbstractState("REFPROP", "MM")
    AS.update(CP.PSmass_INPUTS,p_inlet,entropy)
    T_isobar_inlet = AS.T()
    isobar_inlet.append(T_isobar_inlet)
    AS.update(CP.PSmass_INPUTS,p_outlet,entropy)
    T_isobar_outlet = AS.T()
    isobar_outlet.append(T_isobar_outlet)

# Plot the T-s diagram
plt.figure()
# Plot the colormap of Gamma
contour = plt.contourf(ss, TT, Gamma, 15, cmap='viridis')
# draw a curve through the Gamma_line[i,j] = 1 points
plt.contour(ss, TT, Gamma_line, 1, colors='black', linewidths=0.2)
# Plot the saturated vapor line
plt.plot(s_vaps, Tcrits, 'k',linewidth=0.75)
plt.plot(s_liqs, Tcrits, 'k',linewidth=0.75)
# Plot the isobars
plt.plot(s_gamma, isobar_inlet, 'k--',linewidth=0.75)
plt.plot(s_gamma, isobar_outlet, 'k--',linewidth=0.75)

#add scatter marker to the inlet and outlet of the nozzle
plt.scatter(s[0], T[0], color='red', marker='x')
plt.scatter(s_nozzle_out, T_nozzle_out, color='red', marker='x')
#add text to the inlet (1) and outlet (2) of the nozzle
plt.text(s[0]-10, T[0], '1', ha='right')
plt.text(s_nozzle_out-10, T_nozzle_out, '2', ha='right')
#draw a dashed line between the inlet and outlet of the nozzle
plt.plot([s[0], s_nozzle_out], [T[0], T_nozzle_out], 'k--',linewidth=0.5,color='red')

#plt.title('Colormap of constant Gamma in the T-s diagram')
plt.xlabel('Specific entropy (J/(kg.K))')
plt.xlim(400, s_max)
plt.ylim(T_min, T_max)
plt.ylabel('Temperature (K)')
#plot colorbar named Gamma, relates to contourf, goes from 0.3 to 1.0
cbar = plt.colorbar(contour)
cbar.set_label('Gamma [-]')
#save the plot
plt.savefig('Gamma_T_s_diagram.png', dpi = 600)
plt.show()
"""

#p-v diagram

# Define the pressure and specific volume ranges for the diagram
P_min, P_max = 1e4, 4e6  # Pressure in Pa, modify as needed
v_min, v_max = 1.5e-3, 5e-1  # Specific volume in m^3/kg, modify as needed
steps = 1000  # Number of steps in the grid

# Set up a grid for pressure and specific volume
P_gamma = np.linspace(P_min, P_max, steps)
v_gamma = np.linspace(v_min, v_max, steps)

# Create a meshgrid for P and v
PP, vv = np.meshgrid(P_gamma, v_gamma)

# Initialize an array to hold the state information (e.g., phase)
Gamma = np.zeros_like(PP)
Gamma_line = np.zeros_like(PP)

# Determine Gamma at each point in the grid
for i in range(PP.shape[0]):
    for j in range(PP.shape[1]):
        AS = CP.AbstractState("REFPROP", "MM")  # Use your desired fluid
        #if the next line raises an exception, write NaN
        try:
            AS.update(CP.DmassP_INPUTS, 1/vv[i, j], PP[i, j])
        except:
            Gamma[i, j] = np.nan
            continue
        if AS.Q() >= 0 and AS.Q() <= 1:
            Gamma[i, j] = np.nan  # Assign NaN in two-phase region
        else:
            # Use the CoolProp function to calculate Gamma
            Gamma[i, j] = AS.fundamental_derivative_of_gas_dynamics()
            if Gamma[i, j] > 1:
                Gamma[i, j] = np.nan
#read through the meshgrid and find first points where either NaN switches to a number or vice versa and switch that nan to 1
for i in range(Gamma.shape[0]):
    for j in range(Gamma.shape[1]):
        if np.isnan(Gamma[i, j]) and not np.isnan(Gamma[i-1, j]):
            Gamma_line[i, j] = 1
        if not np.isnan(Gamma[i, j]) and np.isnan(Gamma[i-1, j]):
            Gamma_line[i, j] = 1

#calculate saturated vapor line
v_vaps = []
v_liqs = []
pcrits = []
for Ps in P_gamma:
    if Ps < 1.93e6:
        pcrits.append(Ps)
        AS = CP.AbstractState("REFPROP", "MM")
        AS.update(CP.PQ_INPUTS, Ps,1)
        v_vapor = 1/AS.rhomass() # J/(kg*K)
        AS.update(CP.PQ_INPUTS,Ps,0)
        v_liquid = 1/AS.rhomass()
        v_vaps.append(v_vapor)
        v_liqs.append(v_liquid)

contour = plt.contourf(vv, PP, Gamma, 15, cmap='viridis')
cbar = plt.colorbar(contour)
#add name to the colorbar
cbar.set_label('Gamma [-]')
# draw a curve through the Gamma_line[i,j] = 1 points
#plt.contour(vv, PP, Gamma_line, 1, colors='black', linewidths=0.2)

# Plot the saturated lines
plt.plot(v_vaps, pcrits, 'k',linewidth=0.75)
plt.plot(v_liqs, pcrits, 'k',linewidth=0.75)
# Plot the isobars
p_inlets = []
p_outlets = []
for v in v_gamma:
    p_inlets.append(p_inlet)
    p_outlets.append(p_outlet)
plt.plot(v_gamma, p_inlets, 'k--',linewidth=0.75)
plt.plot(v_gamma, p_outlets, 'k--',linewidth=0.75)
"""
#add scatter marker to the inlet and outlet of the nozzle
plt.scatter(1/rho[0], p[0], color='red', marker='x')
plt.scatter(1/rho_nozzle_out, p_outlet, color='red', marker='x')
#add text to the inlet (1) and outlet (2) of the nozzle
plt.text(1/rho[0]-0.0001, p[0], '1', ha='right')
plt.text(1/rho_nozzle_out-0.0001, p_outlet, '2', ha='right')
#draw a dashed line between the inlet and outlet of the nozzle
plt.plot([1/rho[0], 1/rho_nozzle_out], [p[0], p_outlet], 'k--',linewidth=0.5,color='red')
"""
plt.xlabel(r'Specific Volume (m$^{3}$.kg$^{-1}$)')
plt.ylabel('Pressure (Pa)')
plt.xscale('log')  # Logarithmic scale for better visibility
plt.ylim(P_min, P_max)
plt.xlim(v_min, v_max)
#y axis ticks in scientific notation
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.savefig('Gamma_P_v_diagram.png', dpi = 600)
plt.show()

#compare nozzle design with real gas properties and perfect gas properties - overlap MoC

#calculate the velocity triangle nozzle outlet using alpha_stator

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




