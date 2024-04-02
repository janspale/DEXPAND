from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Load the REFPROP shared library
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])

# Define the REFPROP fluid
fluid = "MM"# hexamethyldisiloxane
z = [1.0]   # Mole fractions
eta_gen = 0.80 # Generator efficiency
eta_vfd = 0.92 # VFD efficiency

#initialize REFPROP
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX']) # Instantiate the REFPROP function library
RP.SETPATHdll(os.environ['RPPREFIX']) # Set the path to the folder containing the REFPROP shared library
RP.SETFLUIDSdll(fluid)  # Set the fluids
baseSI = RP.GETENUMdll(0, "MASSBASESI").iEnum # Get the base SI units
iMass = 0; iFlag = 0 # Set the mass and flag to 0

def check_err(r): # Function to check the error code from REFPROP
    if r.ierr > 0:
        raise ValueError(r.herr)
    
os.chdir('C:\\GitHub\\DEXPAND\\experimental_data\\MM_testrig_CTU\\W50_Spratech') # Change the working directory to the folder containing the data

# Load the experimental data
df = pd.read_csv('w50_data_220324.csv') 


# read values into numpy arrays
time = df['cas'].values

# Convert time from strings in the format 'hh:mm:ss' to time in seconds
df['time_secs'] = df['cas'].apply(lambda t: int(t.split(':')[0])*3600 + int(t.split(':')[1])*60 + int(t.split(':')[2]))

steady_states = []
pressure_levels = [410, 500, 570]  # Define pressure levels
time_ranges = [(10*3600 + 48*60, 11*3600 + 26*60),  # Time ranges for each pressure level
               (13*3600 + 37*60, 14*3600 + 7*60),
               (14*3600 + 36*60, 14*3600 + 55*60)]

for pressure, (start_time, end_time) in zip(pressure_levels, time_ranges):
    # Filter data based on the time range
    filtered_df = df[(df['time_secs'] >= start_time) & (df['time_secs'] <= end_time)]
    
    # Identify unique rpm values within this time range
    unique_rpms = filtered_df['n_gen (rpm)'].unique()
    
    # Initialize a dictionary to hold the average values for each rpm
    rpm_averages = {}
    
    # Calculate averages for each rpm
    for rpm in unique_rpms:
        if 3000 <= rpm <= 8000:  # Considering only rpm values between 3000 and 8000
            rpm_data = filtered_df[filtered_df['n_gen (rpm)'] == rpm]
            averages = {
                'p_ad_avg': rpm_data['p_ad (kPa)'].mean(),
                'p_em_avg': rpm_data['p_kd (kPa)'].mean(),
                'T_ad_avg': rpm_data['t_ad (°C)'].mean(),
                'V_dot_avg': rpm_data['V_MM (l/min)'].mean(),
                'p_nc_avg': rpm_data['p_NC (kPa)'].mean(),
                'T_kd_avg': rpm_data['t_kd (°C)'].mean(),
                'P_gen_avg': rpm_data['P_gen_RM (W)'].mean(),
                'n_gen_avg': rpm_data['n_gen (rpm)'].mean(),
            }
            # Adjust the generator power to get the shaft power
            averages['P_shaft_avg'] = averages['P_gen_avg'] / eta_gen / eta_vfd
            averages['n_shaft_avg'] = averages['n_gen_avg']

            # REFPROP calculations
            T_kd_avg = averages['T_kd_avg'] + 273.15  # Convert to Kelvin
            p_nc_avg = averages['p_nc_avg'] * 1000  # Convert to Pa
            r = RP.REFPROPdll(fluid, "TP", "D", baseSI, iMass, iFlag, T_kd_avg, p_nc_avg, z)
            averages['rho_nc'] = r.Output[0]

            averages['m_dot_avg'] = averages['V_dot_avg'] / 60000 * averages['rho_nc']  # Convert l/min to m^3/s and multiply by density

            # Turbine inlet enthalpy and entropy
            T_ad_avg = averages['T_ad_avg'] + 273.15  # Convert to Kelvin
            p_ad_avg = averages['p_ad_avg'] * 1000  # Convert to Pa
            r = RP.REFPROPdll(fluid, "TP", "H;S", baseSI, iMass, iFlag, T_ad_avg, p_ad_avg, z)
            averages['h_ad_avg'] = r.Output[0]
            averages['s_ad_avg'] = r.Output[1]

            # Isentropic turbine expansion to find outlet enthalpy
            p_em_avg = averages['p_em_avg'] * 1000  # Convert to Pa
            r = RP.REFPROPdll(fluid, "PS", "H", baseSI, iMass, iFlag, p_em_avg, averages['s_ad_avg'], z)
            averages['h_out_is'] = r.Output[0]

            # Turbine efficiency calculation
            averages['eta_turbine'] = averages['P_shaft_avg'] / (averages['m_dot_avg'] * (averages['h_ad_avg'] - averages['h_out_is']))

            rpm_averages[rpm] = averages
    #average the average p_ad_avg 
    #p_ad_avg = np.mean([rpm_averages[rpm]['p_ad_avg'] for rpm in rpm_averages])
    steady_states.append(rpm_averages)
# Plotting
fig, ax = plt.subplots()
ax2 = ax.twinx()
lines, labels = [], []

for i, state in enumerate(steady_states):
    x = [state[rpm]['n_shaft_avg'] for rpm in state]
    y = [state[rpm]['P_shaft_avg'] for rpm in state]
    line1, = ax.plot(x, y, label=f'$P_{{{pressure_levels[i]} kPa}}$')
    lines.append(line1)
    labels.append(line1.get_label())
    ax.scatter(x, y, marker='o')

    # Plotting efficiency on the secondary axis
    w = [state[rpm]['n_shaft_avg'] for rpm in state]
    z = [state[rpm]['eta_turbine'] for rpm in state]
    line2, = ax2.plot(w, z, label=f'$\eta_{{{pressure_levels[i]} kPa}}$', linestyle='--')
    lines.append(line2)
    labels.append(line2.get_label())
    ax2.scatter(w, z, marker='x')

ax2.set_ylabel(r'$\eta_{turbine}$ (-)')
ax2.set_ylim(0, 1)
ax.set_xlabel(r'$n_{shaft}$ (rpm)')
ax.set_ylabel(r'$P_{shaft}$ (W)')
fig.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 1.01), ncol=3)

plt.savefig('P_shaft_vs_n_shaft.png', dpi=600, bbox_inches='tight')
plt.show()

#save the evaluated data to a csv file using pandas
evaluated_data = []
for i, state in enumerate(steady_states):
    for rpm in state:
        evaluated_data.append({
            'pressure': pressure_levels[i],
            'n_shaft': state[rpm]['n_shaft_avg'],
            'P_shaft': state[rpm]['P_shaft_avg'],
            'm_dot': state[rpm]['m_dot_avg'],
            'eta_turbine': state[rpm]['eta_turbine']
        })
evaluated_df = pd.DataFrame(evaluated_data)
evaluated_df.to_csv('evaluated_data_2202324.csv', index=False)

#print the evaluated data to the console
#print(evaluated_df)
