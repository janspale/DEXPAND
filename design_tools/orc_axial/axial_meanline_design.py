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
h_over_D_mid = 0.1  # [-]

#initialize REFPROP
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX']) # Instantiate the REFPROP function library
RP.SETPATHdll(os.environ['RPPREFIX']) # Set the path to the folder containing the REFPROP shared library
RP.SETFLUIDSdll(fluid)  # Set the fluids
baseSI = RP.GETENUMdll(0, "MASSBASESI").iEnum # Get the base SI units
iMass = 0; iFlag = 0 # Set the mass and flag to 0

def check_err(r): # Function to check the error code from REFPROP
    if r.ierr > 0:
        raise ValueError(r.herr)

#initialize state table
h,s,T,p,rho,a,Z,Ma = np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5),np.zeros(5)


# turbine inlet ([0])
p[0] = p_inlet # Pa
T[0] = T_inlet # K
r = RP.REFPROPdll(fluid,"PT","H;S;D;W;Z",baseSI,iMass,iFlag,p_inlet,T_inlet,z); check_err(r)
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
phi_stator = np.sqrt(1-(0.0029*Ma[1]**3-0.0502*Ma[1]**2+0.2241*Ma[1]-0.0877)) 
# [-], velocity loss coefficient at the nozzles

# sonic point (nozzle throat) ([2])


#calculate isentropic expansion quasi-steady 1D, 1000 steps in pressure
