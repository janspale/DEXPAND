# this file is used to evaluate the efficiency of the spindle motor using the experimental data obtained at 08/28/2024 at VUT Brno
# the same spindle is then used in generator mode for the turbine, with the same VFD

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

os.chdir('C:\\GitHub\\DEXPAND\\experimental_data\\MM_testrig_CTU\\VUT_spindle') # Change the working directory to the folder containing the data

# load the data
data = pd.read_csv('py_import_data_spindle.csv')

#read each column into a separate numpy array named after the column
torque = data['Torque [Nm]'].values # Torque in Nm
#Torque [Nm]	f [Hz]	rpm	P_avgd_VFD [W]	P_act_VFD [W]	P_yoko [W]	P_mech [W]	PF [-]	eta_el_VFD [-]	eta_el_yoko [-]	slip [-]	delta_yoko [W]	delta_yoko [-]
f = data['f [Hz]'].values # Frequency in Hz
rpm = data['rpm'].values # Speed in rpm
P_avgd_VFD = data['P_avgd_VFD [W]'].values # Average power in W
P_act_VFD = data['P_act_VFD [W]'].values # Actual power in W
P_yoko = data['P_yoko [W]'].values # Yoko power in W
P_mech = data['P_mech [W]'].values # Mechanical power in W
PF = data['PF [-]'].values # Power factor
eta_el_VFD = data['eta_el_VFD [-]'].values # Electrical efficiency of the VFD
eta_el_yoko = data['eta_el_yoko [-]'].values # Electrical efficiency of the Yoko spindle
slip = data['slip [-]'].values # Slip
delta_yoko = data['delta_yoko [W]'].values # Difference in power between VFD and Yoko spindle
delta_yoko_perc = data['delta_yoko [-]'].values # Difference in power between VFD and Yoko spindle in percentage


# plot the electrical efficiency of the VFD and yoko as a function of torque
plt.figure()
plt.plot(torque,eta_el_VFD, label='VFD')
plt.plot(torque,eta_el_yoko, label='Yoko',linestyle='--')
plt.xlabel('Torque [Nm]')
plt.ylabel('Electrical efficiency [-]')
plt.title('Electrical efficiency as a function of torque')
#polyfit the data
p = np.polyfit(torque,eta_el_VFD,5)
r = np.polyfit(torque,eta_el_yoko,5)
plt.plot(torque,np.polyval(p,torque),'r')
plt.plot(torque,np.polyval(r,torque),'g')
plt.show()

# plot a heatmap of eta_el_VFD as a function of torque and frequency
# and also plot isocontours of eta_el_VFD [0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]

#create a meshgrid of torque and frequency
f_unique = np.unique(f)
torque_unique = np.unique(torque)
f_mesh,torque_mesh = np.meshgrid(f_unique,torque_unique)

# create a meshgrid of eta_el_VFD
eta_el_VFD_mesh = np.zeros((len(torque_unique),len(f_unique)))
for i in range(len(torque_unique)):
    for j in range(len(f_unique)):
        eta_el_VFD_mesh[i,j] = np.mean(eta_el_VFD[(f==f_unique[j]) & (torque==torque_unique[i])])

# plot the heatmap, iso efficiency lines and the values of the iso efficiency
plt.figure()
plt.pcolormesh(f_mesh,torque_mesh,eta_el_VFD_mesh,cmap='viridis')
plt.colorbar()
plt.xlabel('Frequency [Hz]')
plt.ylabel('Torque [Nm]')
plt.title('Electrical efficiency of the VFD')
# also plot iso efficiency lines and the values of the iso efficiency
iso_eff = [0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95]
for eff in iso_eff:
    plt.contour(f_mesh,torque_mesh,eta_el_VFD_mesh,levels=[eff],colors='r',linestyles='--')

# plot the iso efficiency values
for eff in iso_eff:
    idx = np.argmin(np.abs(eta_el_VFD_mesh-eff))
    plt.text(f_mesh.flatten()[idx],torque_mesh.flatten()[idx],str(eff),color='r')

plt.show()


