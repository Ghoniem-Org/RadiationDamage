import pandas as pd
import numpy as np

def read_data():
    # Load C:\Users\Owner\Documents\Repos\RadiationDamage\Creep\material_data.xlsxthe Excel file
    df_steel = pd.read_excel(r'./material_data.xlsx', sheet_name='steel', engine='openpyxl')

    # List of parameter names you want to extract
    input_data = ['tau_oct', 'b', 'D_o', 'SFE', 'delta', 'sigma_o', 'E', 
              'nu', 'N_p', 'r_p', 'k', 'Omega', 'W_g', 'a1', 'c_jog', 
              'E_core', 'E_s', 'E_m', 'K_c', 'Beta', 'xi', 'zeta']

    # Initialize a dictionary to hold parameter values
    p_value = {}

    # Loop through each parameter, extract its value, and store it in the dictionary
    for p in input_data:
        value = df_steel.loc[df_steel['Parameter'] == p, 'Value'].values
        p_value[p] = float(value[0]) if len(value) > 0 else None  # Assign None if the parameter is not found

    # Now you can access each parameter value like so:
    tau_oct = p_value['tau_oct']
    b = p_value['b']
    D_o = p_value['D_o']
    SFE = p_value['SFE']
    delta = p_value['delta']
    sigma_o = p_value['sigma_o']
    E = p_value['E']
    nu = p_value['nu']
    N_p = p_value['N_p']
    r_p = p_value['r_p']
    k = p_value['k']
    Omega = p_value['Omega']
    W_g = p_value['W_g']
    a1 = p_value['a1']
    c_jog = p_value['c_jog']
    E_core = p_value['E_core']
    E_s = p_value['E_s']
    E_m = p_value['E_m']
    K_c = p_value['K_c']
    Beta = p_value['Beta']
    xi = p_value['xi']
    zeta = p_value['zeta']

    # set print precision digits
    np.set_printoptions(precision=10)
    print('tau_oct =', tau_oct)
    print('b =', b)
    print('D_o =', D_o)
    print('SFE =', SFE)
    print('delta =', delta)
    print('sigma_o =', sigma_o)
    print('E =', E)
    print('nu =', nu)
    print('N_p =', N_p)
    print('r_p =', r_p)
    print('k =', k)
    print('Omega =', Omega)
    print('W_g =', W_g)
    print('a1 =', a1)
    print('c_jog =', c_jog)
    print('E_core =', E_core)
    print('E_s =', E_s)
    print('E_m =', E_m)
    print('K_c =', K_c)
    print('Beta =', Beta)
    print('xi =', xi)
    print('zeta =', zeta)
    
    return p_value

def parse_output(output):
    lines = output.strip().split('\n')
    data = [list(map(float, line.strip().split())) for line in lines]
    return np.array(data)
