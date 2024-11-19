"""
This is a Python design tool for the design of an axial impulse single stage small scale ORC turbine.
Accompanies the PhD thesis of J.Spale, 2024. CTU in Prague, supervisor: prof.Ing.M.Kolovratnik, CSc.
for questions, send an email to: Jan.Spale@cvut.cz
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
from deap import base, creator, tools, algorithms
from multiprocessing import Pool
import random

from matplotlib.patches import Arc
from matplotlib.transforms import Bbox, IdentityTransform, TransformedBbox

class AngleAnnotation(Arc):
    """
    Draws an arc between two vectors which appears circular in display space.
    """
    def __init__(self, xy, p1, p2, size=75, unit="points", ax=None,
                 text="", textposition="inside", text_kw=None, **kwargs):
        """
        Parameters
        ----------
        xy, p1, p2 : tuple or array of two floats
            Center position and two points. Angle annotation is drawn between
            the two vectors connecting *p1* and *p2* with *xy*, respectively.
            Units are data coordinates.

        size : float
            Diameter of the angle annotation in units specified by *unit*.

        unit : str
            One of the following strings to specify the unit of *size*:

            * "pixels": pixels
            * "points": points, use points instead of pixels to not have a
              dependence on the DPI
            * "axes width", "axes height": relative units of Axes width, height
            * "axes min", "axes max": minimum or maximum of relative Axes
              width, height

        ax : `matplotlib.axes.Axes`
            The Axes to add the angle annotation to.

        text : str
            The text to mark the angle with.

        textposition : {"inside", "outside", "edge"}
            Whether to show the text in- or outside the arc. "edge" can be used
            for custom positions anchored at the arc's edge.

        text_kw : dict
            Dictionary of arguments passed to the Annotation.

        **kwargs
            Further parameters are passed to `matplotlib.patches.Arc`. Use this
            to specify, color, linewidth etc. of the arc.

        """
        self.ax = ax or plt.gca()
        self._xydata = xy  # in data coordinates
        self.vec1 = p1
        self.vec2 = p2
        self.size = size
        self.unit = unit
        self.textposition = textposition

        super().__init__(self._xydata, size, size, angle=0.0,
                         theta1=self.theta1, theta2=self.theta2, **kwargs)

        self.set_transform(IdentityTransform())
        self.ax.add_patch(self)

        self.kw = dict(ha="center", va="center",
                       xycoords=IdentityTransform(),
                       xytext=(0, 0), textcoords="offset points",
                       annotation_clip=True)
        self.kw.update(text_kw or {})
        self.text = ax.annotate(text, xy=self._center, **self.kw)

    def get_size(self):
        factor = 1.
        if self.unit == "points":
            factor = self.ax.figure.dpi / 72.
        elif self.unit[:4] == "axes":
            b = TransformedBbox(Bbox.unit(), self.ax.transAxes)
            dic = {"max": max(b.width, b.height),
                   "min": min(b.width, b.height),
                   "width": b.width, "height": b.height}
            factor = dic[self.unit[5:]]
        return self.size * factor

    def set_size(self, size):
        self.size = size

    def get_center_in_pixels(self):
        """return center in pixels"""
        return self.ax.transData.transform(self._xydata)

    def set_center(self, xy):
        """set center in data coordinates"""
        self._xydata = xy

    def get_theta(self, vec):
        vec_in_pixels = self.ax.transData.transform(vec) - self._center
        return np.rad2deg(np.arctan2(vec_in_pixels[1], vec_in_pixels[0]))

    def get_theta1(self):
        return self.get_theta(self.vec1)

    def get_theta2(self):
        return self.get_theta(self.vec2)

    def set_theta(self, angle):
        pass

    # Redefine attributes of the Arc to always give values in pixel space
    _center = property(get_center_in_pixels, set_center)
    theta1 = property(get_theta1, set_theta)
    theta2 = property(get_theta2, set_theta)
    width = property(get_size, set_size)
    height = property(get_size, set_size)

    # The following two methods are needed to update the text position.
    def draw(self, renderer):
        self.update_text()
        super().draw(renderer)

    def update_text(self):
        c = self._center
        s = self.get_size()
        angle_span = (self.theta2 - self.theta1) % 360
        angle = np.deg2rad(self.theta1 + angle_span / 2)
        r = s / 2
        if self.textposition == "inside":
            r = s / np.interp(angle_span, [60, 90, 135, 180],
                                          [3.3, 3.5, 3.8, 4])
        self.text.xy = c + r * np.array([np.cos(angle), np.sin(angle)])
        if self.textposition == "outside":
            def R90(a, r, w, h):
                if a < np.arctan(h/2/(r+w/2)):
                    return np.sqrt((r+w/2)**2 + (np.tan(a)*(r+w/2))**2)
                else:
                    c = np.sqrt((w/2)**2+(h/2)**2)
                    T = np.arcsin(c * np.cos(np.pi/2 - a + np.arcsin(h/2/c))/r)
                    xy = r * np.array([np.cos(a + T), np.sin(a + T)])
                    xy += np.array([w/2, h/2])
                    return np.sqrt(np.sum(xy**2))

            def R(a, r, w, h):
                aa = (a % (np.pi/4))*((a % (np.pi/2)) <= np.pi/4) + \
                     (np.pi/4 - (a % (np.pi/4)))*((a % (np.pi/2)) >= np.pi/4)
                return R90(aa, r, *[w, h][::int(np.sign(np.cos(2*a)))])

            bbox = self.text.get_window_extent()
            X = R(angle, r, bbox.width, bbox.height)
            trans = self.ax.figure.dpi_scale_trans.inverted()
            offs = trans.transform(((X-s/2), 0))[0] * 72
            self.text.set_position([offs*np.cos(angle), offs*np.sin(angle)])

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
alpha_stator = 13 # deg
beta_rotor = 23.5  # deg
h_over_D_mid = 0.04  # [-]
eta_guess = 0.7  # [-]
e = 1 # [-] partial admission factor

def meanline_design(D_mid, n_design, alpha_stator, h_over_D_mid, eta_guess, p_inlet, T_inlet, beta_rotor, m_dot):
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

    #calculate isentropic expansion quasi-steady 1D, 1000 steps in pressure
    throat_trigger = False
    p_is = np.linspace(p[0],p[1],1000)
    T_is,rho_is,h_is,a_is,c_is,Ma_is,PR_is,v_is = np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is))
    T_act,rho_act,h_act,a_act,c_act,Ma_act,PR_act,s_act = np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is)),np.zeros(len(p_is))
    for i in range(len(p_is)):
        r = RP.REFPROPdll(fluid,"PS","T;D;H;W",baseSI,iMass,iFlag,p_is[i],s[0],z); check_err(r)
        T_is[i],rho_is[i],h_is[i],a_is[i] = r.Output[0:4] # K, kg/m3, J/kg, m/s
        if h_is[i] < h[0]:
            c_is[i] = np.sqrt(2*(h[0]-h_is[i])) # m/s
        else:
            c_is[i] = 0 # m/s
        Ma_is[i] = c_is[i]/a_is[i] # [-]
        PR_is[i] = p_is[i]/p_outlet # [-]
        v_is[i] = 1/rho_is[i] # m3/kg
        #print(f"Ma_is: {Ma_is[i]}")
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
        if h_act[i] < h[0]:
            c_act[i] = np.sqrt(2*(h[0]-h_act[i])) # m/s
        else:
            c_act[i] = 0 # m/s
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

    #calculate nozzle area, throat, outlet and between
    A_outlet = m_dot/(rho_nozzle_out*c_nozzle_out) # m2
    A_throat = m_dot/(rho_is[throat_index]*c_is[throat_index]) # m2
    A_ratio = A_outlet/A_throat # [-]

    h_over_D_mid = abs(h_over_D_mid) # [-]

    #nozzle design
    ht_nozzle = h_over_D_mid*D_mid # m
    b_nozzles_out = A_outlet/ht_nozzle # m
    b_nozzles_throat = A_throat/ht_nozzle # m
    no_nozzles = 10 # [-] #CONSTANT!!!!! - to be varied
    b_nozzle_out = b_nozzles_out/no_nozzles # m
    b_nozzle_throat = b_nozzles_throat/no_nozzles # m

     #gammas
    """

    CP.set_config_string(CP.ALTERNATIVE_REFPROP_PATH, os.getenv('RPPREFIX'))

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

    #add scatter marker to the inlet and outlet of the nozzle
    plt.scatter(1/rho[0], p[0], color='red', marker='x')
    plt.scatter(1/rho_nozzle_out, p_outlet, color='red', marker='x')
    #add text to the inlet (1) and outlet (2) of the nozzle
    plt.text(1/rho[0]-0.0001, p[0], '1', ha='right')
    plt.text(1/rho_nozzle_out-0.0001, p_outlet, '2', ha='right')
    #draw a dashed line between the inlet and outlet of the nozzle
    plt.plot([1/rho[0], 1/rho_nozzle_out], [p[0], p_outlet], 'k--',linewidth=0.5,color='red')

    plt.xlabel(r'Specific Volume (m$^{3}$.kg$^{-1}$)')
    plt.ylabel('Pressure (Pa)')
    plt.xscale('log')  # Logarithmic scale for better visibility
    plt.ylim(P_min, P_max)
    plt.xlim(v_min, v_max)
    #y axis ticks in scientific notation
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.savefig('Gamma_P_v_diagram.png', dpi = 600)
    plt.show()
    """

    #velocity triangles
    c1a = c_nozzle_out*np.sin(np.radians(alpha_stator)) # m/s; axial absolute velocity stator outlet
    c1u = c_nozzle_out*np.cos(np.radians(alpha_stator)) # m/s; tangential absolute velocity stator outlet
    c1 = np.sqrt(c1a**2+c1u**2) # m/s; absolute velocity stator outlet
    w1a = c1a # m/s; axial relative velocity stator outlet
    w1u = c1u - U # m/s; tangential relative velocity stator outlet
    w1 = np.sqrt(w1a**2+w1u**2) # m/s; relative velocity stator outlet
    Ma1r = w1a/a_nozzle_out # [-]; relative Mach number stator outlet

    beta1 = beta_rotor # deg; rotor inlet flow angle                          
    beta2 = 180 - beta1 # deg; rotor outlet flow angle (impulse turbine) 

    theta = beta2 - beta1 # deg; rotor total deflection angle

    ht_rotor = ht_nozzle + 0.0003 # m; rotor height
    P_guess = m_dot*dh_stage*eta_guess # W; guess for power output
    chords = [0.005,0.0075,0.01,0.0125,0.015,0.2,0.25,0.3,0.4,0.5] #list of chord values
    chord_guess = 1.6*np.sqrt((P_guess/m_dot/1000+25))/1000 # m; guess for rotor chord
    #chord = min(chords, key=lambda x:abs(x-chord_guess)) #pick the closest chord value from the list
    chord = 0.0125 # m; rotor chord
    phi_it = 0.957 - 0.000362 * theta - 0.0258 * Ma1r + 0.00000639 * theta**2 + 0.0674 * Ma1r**2 - 0.0000000753 * theta**3 - 0.043 * Ma1r**3 - 0.000238 * theta * Ma1r + 0.00000145 * theta**2 * Ma1r + 0.0000425 * theta * Ma1r**2 # [-]; impulse turbine loss factor, correlation
    def get_phi(phi_it,ht_rotor,chord):
        #need phi_it loop here to get w2 
        current_phi = phi_it
        k1 = 2 # [-]; phi_it_coeff m
        k2 = 0.65 # [-]; phi_it_coeff n
        phi_corrected = np.sqrt(1-((1-current_phi**2)/k1)*(1+(k1-1)*(ht_rotor/chord)**(-k2))) # [-]; corrected impulse turbine loss factor
        phi_diff = abs(phi_corrected-current_phi) # [-]; difference between corrected and initial impulse turbine loss factor
        return phi_corrected

    phi_rotor = np.real(get_phi(phi_it,ht_rotor,chord)) # [-]; corrected impulse turbine loss factor

    w2 = -w1 * phi_rotor # m/s; relative velocity rotor outlet
    w2a = w2*(-np.sin(np.radians(float(beta2)))) # m/s; axial relative velocity rotor outlet
    w2u = w2*(-np.cos(np.radians(float(beta2)))) # m/s; tangential relative velocity rotor outlet
    #print("w2a:", w2a, "w2u:", w2u)
    if w2.imag != 0:
        print("Warning: complex number for w2")
    c2a = w2a # m/s; axial absolute velocity rotor outlet
    c2u = w2u + U # m/s; tangential absolute velocity rotor outlet
    c2 = np.sqrt(c2a**2+c2u**2) # m/s; absolute velocity rotor outlet
    #print("c2a:", c2a, "c2u:", c2u)
    #if c2 is a complex number, print a warning
    if c2.imag != 0:
        print("Warning: complex number for c2")
    alpha_rotor = np.degrees(np.arctan(float(c2a) / float(c2u))) # deg; rotor outlet absolute angle
    dh_rotor = (w1**2-w2**2)/2 # J/kg; enthalpy drop rotor
    h2 = h_act[-1] + dh_rotor # J/kg; enthalpy rotor outlet
    r = RP.REFPROPdll(fluid,"HP","T;D;W;S",baseSI,iMass,iFlag,h2,p_outlet,z); check_err(r)
    T2,rho2,a2,s2 = r.Output[0:4] # K, kg/m3, m/s, J/kgK
    Ma2 = c2/a2 # [-]; Mach number rotor outlet
    Ma2r = w2/a2 # [-]; relative Mach number rotor outlet
    h2t = h2 + (c2**2)/2 # J/kg; total enthalpy rotor outlet
    #check if the state is superheated vapor
    try:
        q_tot_out = RP.REFPROPdll(fluid,"HS","QMASS",baseSI,iMass,iFlag,h2t,s2,z); check_err(q_tot_out)
        #print(f"q_tot_out: {q_tot_out.Output[0]}")
        if q_tot_out.Output[0] < 0:
            r = RP.REFPROPdll(fluid,"HS","T;D;W;P",baseSI,iMass,iFlag,h2t,s2 ,z); check_err(r)
            T2t,rho2t,a2t,p_outlet_total = r.Output[0:4] # K, kg/m3, m/s, J/kgK
    except:
        T2t,rho2t,a2t,p_outlet_total = T2,rho2,a2,p_outlet # K, kg/m3, m/s, J/kgK
    P_turb_aero = m_dot*U*(c1u-c2u) # W; aerodynamic power
    #print (f"p_outlet_total: {p_outlet_total}")
    #friction loss correlation
    P_fric = 0.01*(n_design/60)**3*D_mid**5*rho2 # W; friction loss
    P_partial_admission = (1-e)*rho2*(n_design/60)**3*D_mid**4*4.5*ht_rotor*1.85/2 # W; partial admission and ventilation loss
    P_mech = P_turb_aero - P_fric - P_partial_admission # W; mechanical power
    eta_turb = (P_mech)/(m_dot*dh_stage) # [-]; turbine efficiency
    h_out_eta = h[0] - eta_turb*dh_stage # J/kg; outlet enthalpy with efficiency
    r=RP.REFPROPdll(fluid,"HP","S",baseSI,iMass,iFlag,h_out_eta,p_outlet,z); check_err(r)
    s_out_eta = r.Output[0] # J/kgK; outlet entropy with efficiency

    #package the results in a dictionary
    res = {"c1u":c1u,"c1a":c1a,"w1u":w1u,"w1a":w1a,"U":U,"c2u":c2u,"c2a":c2a,"w2u":w2u,"w2a":w2a,"alpha_stator":alpha_stator,
           "beta2":beta2,"alpha_rotor":alpha_rotor,"s":s,"h":h,"s_nozzle_out":s_nozzle_out,"h_nozzle_out":h_nozzle_out,"s2":s2,
           "h2":h2,"h2t":h2t,"s_out_eta":s_out_eta,"h_out_eta":h_out_eta,"p_inlet":p_inlet,"p_outlet":p_outlet,
           "p_outlet_total":p_outlet_total,"fluid":fluid,"z":z,"PR_is":PR_is,"Ma_is":Ma_is,"p_is":p_is,"throat_index":throat_index,
           "PR_act":PR_act,"Ma_act":Ma_act,"eta_turb":eta_turb,"P_turb_aero":P_turb_aero,"P_mech":P_mech,"U_over_c_is":U_over_c_is,"Ma1r":Ma1r,
    }

    return res

def draw_velocity_triangles(c1u,c1a,w1u,w1a,U,c2u,c2a,w2u,w2a,alpha_stator,beta2,alpha_rotor):
    plt.figure()
    #====stator outlet====
    #end the line with an arrow
    plt.arrow(0,0,c1u,c1a, head_width=5, head_length=5, fc='red', ec='red', length_includes_head=True)
    plt.arrow(0,0,w1u,w1a, head_width=5, head_length=5, fc='blue', ec='blue', length_includes_head=True)
    plt.arrow(w1u,w1a,U,0, head_width=5, head_length=5, fc='green', ec='green', length_includes_head=True)
    #add text to the arrows
   # plt.text(c1u/2+40,w1a/2+5, r'$\overrightarrow{c_2}$', ha='center', va='center', color='red')
    #plt.text(w1u/2-5,w1a/2+5, r'$\overrightarrow{w_2}$', ha='center', va='center', color='blue')
    #plt.text(w1u+U/2,w1a-4, r'$\overrightarrow{U}$', ha='center', va='center', color='green')
    #draw circular arc for alpha_stator between x axis and c1
    #alpha2 = Arc((0,0),20,10,0,0,alpha_stator,linewidth=1.5, color='red')
    #plt.gca().add_patch(alpha2)
   # plt.text(28,3, r'$\alpha_{2}$', ha='center', va='center', color='red')
    #draw circular arc for beta2 between x axis and w2
    #beta_stator = Arc((0,0),40,20,0,0,180-beta2,linewidth=1.5, color='blue')
    #plt.gca().add_patch(beta_stator)
    #plt.text(30,10, r'$\beta_{2}$', ha='center', va='center', color='blue')
    #====rotor outlet====
    #end the line with an arrow
    plt.arrow(0,0,c2u,c2a, head_width=5, head_length=5, fc='red', ec='red', length_includes_head=True)
    plt.arrow(0,0,w2u,w2a, head_width=5, head_length=5, fc='blue', ec='blue', length_includes_head=True)
    plt.arrow(w2u,w2a,U,0, head_width=5, head_length=5, fc='green', ec='green', length_includes_head=True)
    #add text to the arrows
    #plt.text(c2u/2-9,c2a/2+4, r'$\overrightarrow{c_3}$', ha='center', va='center', color='red')
    #plt.text(w2u/2-25,w2a/2+4, r'$\overrightarrow{w_3}$', ha='center', va='center', color='blue')
   # plt.text(w2u+U/2,w2a-4, r'$\overrightarrow{U}$', ha='center', va='center', color='green')
    #draw circular arc for alpha_rotor between x axis and c2
   # alpha3 = Arc((0,0),20,10,0,0,alpha_rotor,linewidth=1.5, color='red')
   # plt.gca().add_patch(alpha3)
   # plt.text(8,7, r'$\alpha_{3}$', ha='center', va='center', color='red')
    #draw circular arc for beta3 between x axis and w3
   # beta_rotor = Arc((0,0),40,20,0,0,beta2,linewidth=1.5, color='blue')
    #plt.gca().add_patch(beta_rotor)
   # plt.text(-13,11, r'$\beta_{3}$', ha='center', va='center', color='blue')
    plt.ylabel(r'Axial velocity [m.s$^{-1}$]')
    plt.xlabel(r'Tangential velocity [m.s$^{-1}$]')
    plt.xlim(-150,300)
    plt.ylim(0,100)
    #move x axis to y=0
    plt.axvline(0, color='black',linewidth=0.5)
    #invert y axis
    plt.gca().invert_yaxis()
    #show grid
    plt.grid(True, linestyle='--')
    #equal aspect ratio
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig('velocity_triangles.png', dpi = 1000)
    plt.show()
    plt.close()

def draw_expansion_line(s,h,s_nozzle_out,h_nozzle_out,s2,h2,h2t,s_out_eta,h_out_eta,p_inlet,p_outlet,p_outlet_total,fluid,z):
    #initialize REFPROP
    RP = REFPROPFunctionLibrary(os.environ['RPPREFIX']) # Instantiate the REFPROP function library
    RP.SETPATHdll(os.environ['RPPREFIX']) # Set the path to the folder containing the REFPROP shared library
    RP.SETFLUIDSdll(fluid)  # Set the fluids
    baseSI = RP.GETENUMdll(0, "MASSBASESI").iEnum # Get the base SI units
    iMass = 0; iFlag = 0 # Set the mass and flag to 0

    def check_err(r): # Function to check the error code from REFPROP
        if r.ierr > 0:
            raise ValueError(r.herr)
    s_range = np.linspace(s[0]-20,s_nozzle_out+30,1000)
    isobar_inlet = []
    isobar_outlet = []
    isobar_rotor_out_total = []
    for ss in s_range:
        r = RP.REFPROPdll(fluid,"PS","H",baseSI,iMass,iFlag,p_inlet,ss,z); check_err(r)
        hh = r.Output[0] # J/kg
        isobar_inlet.append(hh)
        r = RP.REFPROPdll(fluid,"PS","H",baseSI,iMass,iFlag,p_outlet,ss,z); check_err(r)
        hh = r.Output[0] # J/kg
        isobar_outlet.append(hh)
        r = RP.REFPROPdll(fluid,"PS","H",baseSI,iMass,iFlag,p_outlet_total,ss,z); check_err(r)
        hh = r.Output[0] # J/kg
        isobar_rotor_out_total.append(hh)

    fig, ax = plt.subplots()
    ax.plot(s_range,isobar_inlet, label='Isobar inlet', color='black', linestyle='--', linewidth = 0.7)
    ax.plot(s_range,isobar_outlet, label='Isobar outlet', color='black', linestyle='--', linewidth = 0.7)
    ax.plot(s_range,isobar_rotor_out_total, label='Isobar rotor outlet total', color='black', linestyle='--', linewidth = 0.7)
    ax.scatter(s[0], h[0], color='red', marker='x')
    ax.scatter(s_nozzle_out, h_nozzle_out, color='red', marker='x')
    #connect s0 and snozzle out with a dashed line
    ax.plot([s[0], s_nozzle_out], [h[0], h_nozzle_out],linewidth=1,color='black')
    ax.scatter(s[1], h[1], color='red', marker='x')
    ax.scatter(s2, h2t, color='red', marker='x')
    ax.plot([s_nozzle_out, s2], [h_nozzle_out, h2],linewidth=1,color='black')
    ax.scatter(s2, h2, color='red', marker='x')
    ax.plot([s2,s2], [h2,h2t],linewidth=1,color='black')
    ax.scatter(s_out_eta,h_out_eta, color='red', marker='x')
    ax.plot([s[0],s_out_eta], [h[0],h_out_eta],"k--",linewidth=1,color='black')
    #connect s0 and s1 with a dashed line
    ax.plot([s[0], s[1]], [h[0], h[1]], 'k--',linewidth=0.5,color='black')
    ax.grid(True, linestyle='--')
    ax.set(xlabel=r'Specific entropy $s$ [J.(kg.K)$^{-1}$]', ylabel=r'Specific enthalpy $h$ [J.kg$^{-1}$]')
    #change y axis to scientific notation
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #add text to the isobars
    ax.text(s[0]-10, h[0]-2000, r'$p_{evap}$', ha='right')
    ax.text(s_nozzle_out-10, h_nozzle_out-9000, r'$p_{cond}$', ha='right')
    ax.text(s_nozzle_out-10, h_nozzle_out-2000, r'$p_{cond,tot}$', ha='right')
    #add text to the state points
    ax.text(s[0]-1, h[0]+500, '1', ha='right')
    ax.text(s_nozzle_out-1, h_nozzle_out-1000, '2', ha='right')
    ax.text(s[1]-1, h[1]+500, r'2$_{is}$=3$_{is}$', ha='right')
    ax.text(s2-1, h2-3500, '3', ha='right')
    ax.text(s2-1, h2t, r'3$_{t}$', ha='right')
    ax.text(s_out_eta+1, h_out_eta-2000, r'2$_{\eta}$=3$_{\eta}$', ha='left')
    plt.savefig('expansion_line_h_s.png', dpi = 600)
    plt.show()
    plt.close()

def draw_nozzle_expansion(PR_is,Ma_is,p_is,throat_index,PR_act,Ma_act):
    #draw isentropic expansion
    fig, ax = plt.subplots()
    ax.plot(PR_is,Ma_is, label='Isentropic expansion', color='black')
    ax.set(xlabel=r'Pressure ratio $\Pi$ [-]', ylabel=r'Mach number $Ma$ [-]')
    ax.scatter(p_is[throat_index]/p_outlet,1, color='red', marker='x')
    #add text to the throat
    ax.text(p_is[throat_index]/p_outlet+0.05, 1.05, 'Throat', ha='right')
    ax.plot(PR_act,Ma_act, label='Real expansion with losses', color='black', linestyle='--')
    #add scatter to the nozzle outlet for both and add text with the Ma value, use x markers in red
    ax.scatter(p_is[-1]/p_outlet,Ma_act[-1], color='red', marker='x')
    ax.text(p_is[-1]/p_outlet, Ma_act[-1], r'$Ma_2$ = '+str(round(Ma_act[-1],3)), ha='right')
    ax.scatter(p_is[-1]/p_outlet,Ma_is[-1], color='red', marker='x')
    ax.text(p_is[-1]/p_outlet, Ma_is[-1], r'$Ma_{2,is}$ = '+str(round(Ma_is[-1],3)), ha='right')
    #add a tick at x axis for the inlet, throat and for the outlet of the nozzle
    ax.axvline(p_is[0]/p_outlet, color='black', linestyle='--', linewidth=0.5)
    ax.axvline(p_is[throat_index]/p_outlet, color='black', linestyle='--', linewidth=0.5)
    ax.axvline(p_is[-1]/p_outlet, color='black', linestyle='--', linewidth=0.5)
    #add labels to the ticks, rotate the text by 90°, offset the text by 0.5 in x axis
    ax.text(p_is[0]/p_outlet+0.05, 0.5, 'Inlet', ha='right', rotation=90)
    ax.text(p_is[throat_index]/p_outlet+0.05, 0.5, 'Throat', ha='right', rotation=90)
    ax.text(p_is[-1]/p_outlet+0.05, 0.5, 'Outlet', ha='right', rotation=90)
    #write the x axis value next to the tick text, rotated by 90°
    ax.text(p_is[0]/p_outlet-0.05, 0.5, r"$\Pi$ = "+str(round(p_is[0]/p_outlet,2)), ha='left', rotation=90)
    ax.text(p_is[throat_index]/p_outlet-0.05, 0.5, r"$\Pi$ = "+str(round(p_is[throat_index]/p_outlet,2)), ha='left', rotation=90)
    ax.text(p_is[-1]/p_outlet-0.05, 0.5, r"$\Pi$= "+str(round(p_is[-1]/p_outlet,2)), ha='left', rotation=90)
    ax.legend()
    ax.grid(True, linestyle='--')
    #invert x axis
    ax.invert_xaxis()
    plt.savefig('isentropic_expansion_nozzle_losses.png', dpi = 600)
    plt.show()
    plt.close()

def legacy_sensitivity_analysis_rpm(n_design,D_mid,alpha_stator,h_over_D_mid,eta_guess,p_inlet,T_inlet):
    #vary rotational speed in the range of 50 to 150% of nominal speed and calculate the eta_turb, plot it over U_over_c_is - legacy code
    
    n_range = np.linspace(0.5*n_design,1.5*n_design,100)
    eta_turb_range = []
    U_over_c_is_range = []
    for n in n_range:
        res=meanline_design(D_mid,n,alpha_stator,h_over_D_mid,eta_guess,p_inlet,T_inlet)
        eta_turb_range.append(res["eta_turb"])
        U_over_c_is_range.append(res["U_over_c_is"])
    #find maximum eta_turb and corresponding U_over_c_is and rpm
    max_eta_turb = max(eta_turb_range)
    max_eta_turb_index = eta_turb_range.index(max_eta_turb)
    max_eta_turb_rpm = n_range[max_eta_turb_index]
    print(f"Maximum efficiency: {max_eta_turb} at U/c_is: {U_over_c_is_range[max_eta_turb_index]} and rpm: {max_eta_turb_rpm}")
    fig, ax = plt.subplots()
    ax.plot(U_over_c_is_range,eta_turb_range, label='Turbine efficiency', color='black')
    ax.set(xlabel='U/c_is [-]', ylabel='Turbine efficiency [-]')
    #ax.legend()
    ax.grid(True, linestyle='--')
    plt.savefig('sensitivity_analysis_rpm.png', dpi = 600)
    plt.show()
    plt.close()
    
def calculate_efficiency_for_n(args):
    n, D_mid, alpha_stator, h_over_D_mid, eta_guess, p_inlet, T_inlet, m_dot = args
    res = meanline_design(D_mid, n, alpha_stator, h_over_D_mid, eta_guess, p_inlet, T_inlet, beta_rotor, m_dot)
    return res["eta_turb"], res["P_turb_aero"]

def sensitivity_analysis_rpm_pool(n_design, D_mid, alpha_stator, h_over_D_mid, eta_guess, p_inlet, T_inlet):
    p_inlet_range = np.linspace(100000, 700000, 100)  # Pa
    n_range = np.linspace(3000, 20000, 100)           # rpm

    eta_turb_range = []
    P_turb_range = []

    dT_SH = 10  # K; superheating temperature

    # Initialize REFPROP (done outside the loop to avoid reinitialization in each process)
    RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
    RP.SETPATHdll(os.environ['RPPREFIX'])
    RP.SETFLUIDSdll(fluid)
    baseSI = RP.GETENUMdll(0, "MASSBASESI").iEnum
    iMass = 0
    iFlag = 0

    n_cpus = max(1, os.cpu_count() - 4)

    def check_err(r):
        if r.ierr > 0:
            raise ValueError(r.herr)

    for p_inlet in p_inlet_range:
        # Convert p_inlet from Pa to kPa
        p_inlet_kPa = p_inlet / 1000  # Convert Pa to kPa

        # Calculate V̇ in l/min using the given relationship
        V_dot_lpm = (p_inlet_kPa - 14.48) / 29.136  # l/min

        # Check for negative or zero volumetric flow rate
        if V_dot_lpm <= 0:
            continue  # Skip this iteration if V̇ is not physically meaningful

        # Convert V̇ from l/min to m³/s
        V_dot_m3s = V_dot_lpm * (1e-3) / 60  # m³/s

        # Calculate the inlet temperature based on the saturation temperature plus superheat
        r = RP.REFPROPdll(fluid, "PQ", "T", baseSI, iMass, iFlag, p_inlet, 1, z)
        check_err(r)
        T_inlet = r.Output[0] + dT_SH  # K; inlet temperature

        # Calculate the density at the new inlet conditions
        r = RP.REFPROPdll(fluid, "PT", "D", baseSI, iMass, iFlag, p_inlet, 335, z)
        check_err(r)
        rho_inlet = r.Output[0]  # kg/m³

        # Update the mass flow rate based on the new density
        m_dot = V_dot_m3s * rho_inlet  # kg/s

        # Prepare arguments for parallel processing
        args = [(n, D_mid, alpha_stator, h_over_D_mid, eta_guess, p_inlet, T_inlet, m_dot) for n in n_range]

        # Parallel computation
        with Pool(processes=n_cpus) as pool:
            results = pool.map(calculate_efficiency_for_n, args)

        eta_turb_for_p, P_turb_for_p = zip(*results)
        eta_turb_range.append(eta_turb_for_p)
        P_turb_range.append(P_turb_for_p)


    # Plotting the results
    fig, ax = plt.subplots()
    for i, p_inlet in enumerate(p_inlet_range):
        ax.plot(n_range, eta_turb_range[i], label=f'p_inlet = {p_inlet:.0f} Pa')

    ax.set(xlabel='Rotational speed [rpm]', ylabel='Turbine efficiency [-]')
    ax.grid(True, linestyle='--')
    plt.savefig('sensitivity_analysis_rpm.png', dpi=600)
    plt.show()
    plt.close()

    fig, ax = plt.subplots()
    for i, p_inlet in enumerate(p_inlet_range):
        ax.plot(n_range, P_turb_range[i], label=f'p_inlet = {p_inlet:.0f} Pa')
    
    ax.set(xlabel='Rotational speed [rpm]', ylabel='Turbine power [W]')
    ax.grid(True, linestyle='--')
    plt.savefig('sensitivity_analysis_rpm_power.png', dpi=600)
    plt.show()
    plt.close()

    # Find maxima for each P_turb_range and plot it over n
    max_P_turb = [max(P_turb) for P_turb in P_turb_range]
    max_P_turb_index = [P_turb.index(max_P) for P_turb, max_P in zip(P_turb_range, max_P_turb)]
    max_P_turb_rpm = [n_range[index] for index in max_P_turb_index]
    
    # Polyfit the max_P_turb_rpm over p_inlet_range/p_outlet
    p = np.polyfit(p_inlet_range / p_outlet, max_P_turb_rpm, 4)
    
    # Plot the polyfit into the same graph
    ax.plot(p_inlet_range / p_outlet, np.polyval(p, p_inlet_range / p_outlet), 'r--')
    plt.savefig('sensitivity_analysis_rpm_power_max.png', dpi=600)
    plt.show()
    plt.close()

    # Write the polyfit coefficients to the terminal
    print(f"Polyfit coefficients: {p}")

    #create a heatmap of the turbine power over the rotational speed (x) and the pressure ratio (y), use matplotlib
    fig, ax = plt.subplots()
    n_range, p_inlet_range = np.meshgrid(n_range, p_inlet_range)
    P_turb_range = np.array(P_turb_range)
    c = ax.pcolormesh(n_range, p_inlet_range/p_outlet, P_turb_range, cmap='viridis')
    #plot iso-power lines
    iso_power = [1000,2000,3000,4000,5000,6000,7000,8000,9000]
    for power in iso_power:
        ax.contour(n_range, p_inlet_range/p_outlet, P_turb_range, levels=[power], colors='black')
    fig.colorbar(c, ax=ax, label=r'Turbine power $P_{mech}$ [W]')
    ax.set(xlabel=r'Rotational speed $n$ [rpm]', ylabel=r'Pressure ratio $\Pi$ [-]')
    plt.savefig('sensitivity_analysis_rpm_power_heatmap.png', dpi=600)
    plt.show()
    plt.close()

    #create a heatmap of the turbine efficiency over the rotational speed (x) and the pressure ratio (y), use matplotlib, show iso-efficiency lines 0.8,0.775,0.75,0.725,0.7,0.675,0.65,0.625,0.6
    fig, ax = plt.subplots()
    c = ax.pcolormesh(n_range, p_inlet_range/p_outlet, eta_turb_range, cmap='viridis')
    fig.colorbar(c, ax=ax, label=r'Turbine efficiency $\eta_{{is}_{t-s}}$ [-]')
    ax.set(xlabel=r'Rotational speed $n$ [rpm]', ylabel=r'Pressure ratio $\Pi$ [-]')
    #plot iso-efficiency lines
    iso_efficiency = [0.75,0.725,0.7,0.675,0.65,0.625,0.6]
    for eff in iso_efficiency:
        ax.contour(n_range, p_inlet_range/p_outlet, eta_turb_range, levels=[eff], colors='black')
    plt.savefig('sensitivity_analysis_rpm_efficiency_heatmap.png', dpi=600)
    plt.show()
    plt.close()

    #create a heatmap of the turbine torque over the rotational speed (x) and the pressure ratio (y), use matplotlib, show iso-torque lines 2,3,4,5,6,7,8,9,10,11,12 Nm
    fig, ax = plt.subplots()
    c = ax.pcolormesh(n_range, p_inlet_range/p_outlet, P_turb_range/(2*np.pi*n_range/60), cmap='viridis')
    fig.colorbar(c, ax=ax, label=r'Turbine torque $T_{mech}$ [Nm]')
    ax.set(xlabel=r'Rotational speed $n$ [rpm]', ylabel=r'Pressure ratio $\Pi$ [-]')
    #plot iso-torque lines
    iso_torque = [2,3,4,5,6,7,8,9,10,11,12]
    for torque in iso_torque:
        ax.contour(n_range, p_inlet_range/p_outlet, P_turb_range/(2*np.pi*n_range/60), levels=[torque], colors='black')
    plt.savefig('sensitivity_analysis_rpm_torque_heatmap.png', dpi=600)
    plt.show()
    plt.close()

def sensitivity_analysis_alpha_stator(n_design,D_mid,alpha_stator,h_over_D_mid,eta_guess):
    #vary alpha_stator in the range of 10° to 20° and calculate the eta_turb, plot it over alpha_stator
    alpha_stator_range = np.linspace(10,20,100)
    eta_turb_range = []
    for alpha_stator in alpha_stator_range:
        res=meanline_design(D_mid,n_design,alpha_stator,h_over_D_mid,eta_guess)
        eta_turb_range.append(res["eta_turb"])
    fig, ax = plt.subplots()
    ax.plot(alpha_stator_range,eta_turb_range, label='Turbine efficiency', color='black')
    ax.set(xlabel='Stator inlet angle [°]', ylabel='Turbine efficiency [-]')
    #ax.legend()
    ax.grid(True, linestyle='--')
    plt.savefig('sensitivity_analysis_alpha_stator.png', dpi = 600)
    plt.show()
    plt.close()

def sensitivity_analysis_D_mid(n_design,D_mid,alpha_stator,h_over_D_mid,eta_guess):
    #vary alpha_stator in the range of 0.1m to 0.18m and calculate the eta_turb, plot it over D_mid
    D_mid_range = np.linspace(0.1,0.18,100)
    eta_turb_range = []
    for D_mid in D_mid_range:
        res=meanline_design(D_mid,n_design,alpha_stator,h_over_D_mid,eta_guess)
        eta_turb_range.append(res["eta_turb"])
    fig, ax = plt.subplots()
    ax.plot(D_mid_range,eta_turb_range, label='Turbine efficiency', color='black')
    ax.set(xlabel='Mean diameter [m]', ylabel='Turbine efficiency [-]')
    #ax.legend()
    ax.grid(True, linestyle='--')
    plt.savefig('sensitivity_analysis_D_mid.png', dpi = 600)
    plt.show()
    plt.close()

def evaluate(individual):
    D_mid, n_design, alpha_stator, h_over_D_mid,beta_rotor = individual
    result = meanline_design(D_mid, n_design, alpha_stator, h_over_D_mid, eta_guess,p_inlet,T_inlet,beta_rotor, m_dot)
    return (result["eta_turb"],)

def run_GA(population, toolbox, ngen, stats):
    result, log = algorithms.eaSimple(population, toolbox, cxpb=0.5, mutpb=0.2, ngen=ngen, stats=stats, verbose=True)
    return result, log

def plot_convergence(log):
    gen = log.select("gen")
    max_fitness = log.select("max")
    avg_fitness = log.select("avg")

    fig, ax1 = plt.subplots()
    line1 = ax1.plot(gen, max_fitness, "b-", label="Maximum Fitness")
    ax1.set_xlabel("Generation")
    ax1.set_ylabel("Fitness", color="b")
    for tl in ax1.get_yticklabels():
        tl.set_color("b")

    ax2 = ax1.twinx()
    line2 = ax2.plot(gen, avg_fitness, "r-", label="Average Fitness")
    ax2.set_ylabel("Average Fitness", color="r")
    for tl in ax2.get_yticklabels():
        tl.set_color("r")

    lns = line1 + line2
    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc="center right")

    plt.show()

def plot_parameters_evolution(log, param_labels, population):
    gen = log.select("gen")
    num_gens = len(gen)
    fig, axs = plt.subplots(len(param_labels))

    for i, param_label in enumerate(param_labels):
        # Collect the parameter values over the generations
        values_over_gens = [[] for _ in range(num_gens)]
        for g in range(num_gens):
            # Get individuals from the population of the specific generation
            for ind in population:
                values_over_gens[g].append(ind[i])
        # Calculate average parameter values for each generation
        avg_values = [np.mean(values) for values in values_over_gens]
        axs[i].plot(gen, avg_values, label=f'Parameter {param_label}')
        axs[i].set_xlabel("Generation")
        axs[i].set_ylabel(f'{param_label}')
        axs[i].legend()

    plt.tight_layout()
    plt.show()

def plot_fitness_landscape():
    # This requires generating a grid of parameter values and evaluating their fitness
    # Simplified here for illustration
    from mpl_toolkits.mplot3d import Axes3D

    x = np.linspace(0.12, 0.15, 100)
    y = np.linspace(15000, 20000, 100)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(np.sqrt(X**2 + Y**2))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis')

    ax.set_xlabel('D_mid')
    ax.set_ylabel('n_design')
    ax.set_zlabel('Fitness')
    plt.show()

def plot_population_distribution(population):
    D_mid = [ind[0] for ind in population]
    n_design = [ind[1] for ind in population]
    alpha_stator = [ind[2] for ind in population]
    h_over_D_mid = [ind[3] for ind in population]

    fig, axs = plt.subplots(2, 2)
    axs[0, 0].hist(D_mid, bins=20)
    axs[0, 0].set_title('D_mid')
    axs[0, 1].hist(n_design, bins=20)
    axs[0, 1].set_title('n_design')
    axs[1, 0].hist(alpha_stator, bins=20)
    axs[1, 0].set_title('alpha_stator')
    axs[1, 1].hist(h_over_D_mid, bins=20)
    axs[1, 1].set_title('h_over_D_mid')
    fig.tight_layout()
    plt.show()

if __name__=='__main__':
    #run 1DTDT
    tic = timeit.default_timer()

    #senstivity analyses
    #legacy_sensitivity_analysis_rpm(n_design,D_mid,alpha_stator,h_over_D_mid,eta_guess,p_inlet,T_inlet)
    #sensitivity_analysis_rpm_pool(n_design,D_mid,alpha_stator,h_over_D_mid,eta_guess,p_inlet,T_inlet)
    #sensitivity_analysis_alpha_stator(n_design,D_mid,alpha_stator,h_over_D_mid,eta_guess)
    #sensitivity_analysis_D_mid(n_design,D_mid,alpha_stator,h_over_D_mid,eta_guess)
     
    res=meanline_design(D_mid,n_design,alpha_stator,h_over_D_mid,eta_guess,p_inlet,T_inlet,beta_rotor, m_dot)
    print(res)

    #uncomment as you wish to draw the plots
    #draw_nozzle_expansion(res["PR_is"],res["Ma_is"],res["p_is"],res["throat_index"],res["PR_act"],res["Ma_act"])
    #draw_expansion_line(res["s"],res["h"],res["s_nozzle_out"],res["h_nozzle_out"],res["s2"],res["h2"],res["h2t"],res["s_out_eta"],res["h_out_eta"],res["p_inlet"],res["p_outlet"],res["p_outlet_total"],res["fluid"],res["z"])
    draw_velocity_triangles(res["c1u"],res["c1a"],res["w1u"],res["w1a"],res["U"],res["c2u"],res["c2a"],res["w2u"],res["w2a"],res["alpha_stator"],res["beta2"],res["alpha_rotor"])
    #save the variables shown in the velocity triangles into a dictionary
    #velocity_triangles = {"c1u":res["c1u"],"c1a":res["c1a"],"w1u":res["w1u"],"w1a":res["w1a"],"U":res["U"],"c2u":res["c2u"],"c2a":res["c2a"],"w2u":res["w2u"],"w2a":res["w2a"],"alpha_stator":res["alpha_stator"],"beta2":res["beta2"],"alpha_rotor":res["alpha_rotor"]}
    """
    #Genetic Algorithm 
    import multiprocessing
    pool = multiprocessing.Pool(multiprocessing.cpu_count()-2)

    # Define the optimization problem to maximize eta_turb
    creator.create("FitnessMax", base.Fitness, weights=(1.0,)) # Maximization problem
    creator.create("Individual", list, fitness=creator.FitnessMax) # Individual class

    # Setup the parameter ranges
    parameter_bounds = [
        (0.08, 0.2),    # D_mid
        (10000, 25000),  # n_design
        (10, 18),         # alpha_stator (degrees)
        (0.03, 0.07),     # h_over_D_mid
        (15, 30)            #beta_rotor (degrees)
    ]

    param_labels = ['D_mid', 'n_design', 'alpha_stator', 'h_over_D_mid', 'beta_rotor']

    # Helper to create a random individual
    def create_individual():
        return [random.uniform(low, high) for low, high in parameter_bounds]

    # Register the genetic algorithm functions
    toolbox = base.Toolbox()

    toolbox.register("map", pool.map) # Parallel map
    toolbox.register("individual", tools.initIterate, creator.Individual, create_individual) # Individual
    toolbox.register("population", tools.initRepeat, list, toolbox.individual) # Population
    toolbox.register("evaluate", evaluate) # Evaluation function
    toolbox.register("mate", tools.cxBlend, alpha=0.5) # Blend crossover
    toolbox.register("mutate", tools.mutPolynomialBounded, low=[lb for lb, ub in parameter_bounds], up=[ub for lb, ub in parameter_bounds], eta=0.5, indpb=0.1) # Polynomial mutation
    toolbox.register("select", tools.selTournament, tournsize=3) # Tournament selection
    
    # Generate the initial population
    population = toolbox.population(n=40)

    # Run the genetic algorithm
    ngen = 50  # Number of generations
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    print("Population size:", len(population))
    print("Example individual:", population[0])
    print("Crossover probability (cxpb):", 0.5)
    print("Mutation probability (mutpb):", 0.2)
    print("Number of generations (ngen):", ngen)

    result, log = run_GA(population, toolbox, ngen, stats)

    # Find the best result
    best_ind = tools.selBest(result, k=1)[0]
    print('Best Individual: ', best_ind)
    print('Best Fitness: ', best_ind.fitness.values[0])

    pool.close()
    pool.join()

    plot_convergence(log) # Plot the convergence of the GA
    # plot the evolution of different parameters over generations
    #plot_parameters_evolution(log, param_labels,population)
    """
    toc = timeit.default_timer()
    print(f"Elapsed time: {toc-tic} s")