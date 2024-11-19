"""
this code uses REFPROP to draw an h-s diagram of an expansion in a turbine with a given inlet state and a given outlet pressure

Jan Spale, 03/14/2024
Purdue University
jan.spale@cvut.cz
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import os
# pip-installable packages
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
from units import convert

# Set the REFPROP path
RP_path = "c:\\Program Files (x86)\\REFPROP"

# Set the fluids
fluid = "MM" # Set the fluid to be used
zz = [1.0] # Set the mole fractions

#initialize REFPROP
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX']) # Instantiate the REFPROP function library
RP.SETPATHdll(os.environ['RPPREFIX']) # Set the path to the folder containing the REFPROP shared library
RP.SETFLUIDSdll(fluid)  # Set the fluids
baseSI = RP.GETENUMdll(0, "MASSBASESI").iEnum # Get the base SI units
iMass = 0; iFlag = 0 # Set the mass and flag to 0

def check_err(r): # Function to check the error code from REFPROP
    if r.ierr > 0:
        raise ValueError(r.herr)

# Set the inlet state
p_inlet = convert(550,"kPa","Pa") # kPa
T_inlet = convert(180,"C","K") # C
p_outlet = convert(55,"kPa","Pa") # kPa
eta_is = 0.7  # isentropic efficiency

c_in = 10 # m/s, inlet velocity

#calculate inlet state enthalpy
r = RP.REFPROPdll(fluid,"PT","H;S",baseSI,iMass,iFlag,p_inlet,T_inlet,zz); check_err(r)
h_inlet, s_inlet = r.Output[0], r.Output[1]

#calculate outlet state enthalpy isentropic
r = RP.REFPROPdll(fluid,"PS","H;T",baseSI,iMass,iFlag,p_outlet,s_inlet,zz); check_err(r)
h_outlet_is, T_outlet_is = r.Output[0], r.Output[1]

#calculate outlet state enthalpy with real eta
h_outlet = h_inlet - eta_is*(h_inlet-h_outlet_is)
r = RP.REFPROPdll(fluid,"PH","T;S",baseSI,iMass,iFlag,p_outlet,h_outlet,zz); check_err(r)
T_outlet, s_outlet = r.Output[0], r.Output[1]

print("Inlet state: p = {:.2f} kPa, T = {:.2f} C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg/K".format(convert(p_inlet,"Pa","kPa"),convert(T_inlet,"K","C"),h_inlet/1000,s_inlet/1000))
print("Outlet state isentropic: p = {:.2f} kPa, T = {:.2f} C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg/K".format(convert(p_outlet,"Pa","kPa"),T_outlet_is-273.15,h_outlet_is/1000,s_inlet/1000))
print("Outlet state real: p = {:.2f} kPa, T = {:.2f} C, h = {:.2f} kJ/kg, s = {:.2f} kJ/kg/K".format(convert(p_outlet,"Pa","kPa"),T_outlet-273.15,h_outlet/1000,s_outlet/1000))