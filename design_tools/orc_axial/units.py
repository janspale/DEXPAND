"""
perform units conversion

@ Jan Spale, Purdue University, 2023
"""
"""
x = convert(x, a, b)

Inputs:
    x: parameter value
    a: original unit
    b: new unit
    
Outputs:
    x: parameter value after unit conversion 
"""


def convert(x, a, b):
    if a == 'C' and b == 'K':
        x += 273.15
    elif a == 'K' and b == 'C':
        x -= 273.15
    elif a == 'C' and b == 'F':
        x = (x*9/5)+32
    elif a == 'F' and b == 'C':
        x = (x-32)*5/9
    elif a == 'K' and b == 'F':
        x = (x-273.15)*9/5+32
    elif a == 'F' and b == 'K':
        x = (x-32)*5/9+273.15
    elif a == 'delta_K' and b == 'delta_F':
        x = x*1.8
    elif a == 'delta_F' and b == 'delta_K':
        x = x*5/9   
    elif (a == 'kW' and b == 'W') or (a == 'kPa' and b == 'Pa') or (a == 'kJ/kg' and b == 'J/kg') or (a == 'kJ/kg-K' and b == 'J/kg-K'):
        x = x*1000
    elif (a == 'W' and b == 'kW') or (a == 'Pa' and b == 'kPa') or (a == 'J/kg' and b == 'kJ/kg') or (a == 'J/kg-K' and b == 'kJ/kg-K'):
        x = x/1000
    elif a == 'kg/hr' and b == 'kg/s':
        x = x/3600
    elif a == 'kg/s' and b == 'kg/hr':
        x = x*3600
    elif a == 'kW' and b == 'Btu/hr':
        x = x*3412.14245
    elif a == 'Btu/hr' and b == 'kW':
        x = x/3412.14245
    elif a == 'W' and b == 'Btu/hr':
        x = x*3.41214245
    elif a == 'Btu/hr' and b == 'W':
        x = x/3.41214245
    elif a == 'rpm' and b == 'rps':
        x = x/60
    elif a == 'in3/rev' and b == 'm3/rev':
        x = x*1.6387064e-5
    elif a == 'cfm' and b == 'm3/s':
        x = x*4.71947e-4
    elif a == 'm3/s' and b == 'cfm':
        x = x/4.71947e-4
    else:
        raise ValueError(f'unit conversion is undefined')

    return x