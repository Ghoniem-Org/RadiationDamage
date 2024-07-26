import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import subprocess


def read_data():
    # Load the Excel file
    df_steel = pd.read_excel('material_data.xlsx', sheet_name='steel', engine='openpyxl')

    # List of parameter names you want to extract
    input_data = ['nu_v', 'nu_i', 'nu_g', 'Em_v', 'Em_i', 'Em_g',
                  'Eb_v_g', 'Eb_v_2g', 'Eb_2g', 'Ef_v', 'a0', 'Omega',
                  'f', 'b', 'k_B', 'gamma_b', 'B', 'r_ppt', 'd',
                  'N_ppt', 'rho', 'Z_i']

    # Initialize a dictionary to hold parameter values
    p_value = {}

    # Loop through each parameter, extract its value, and store it in the dictionary
    for p in input_data:
        value = df_steel.loc[df_steel['Parameter'] == p, 'Value'].values
        p_value[p] = float(value[0]) if len(value) > 0 else None  # Assign None if the parameter is not found

    # Now you can access each parameter value like so:
    nu_v = p_value['nu_v']
    nu_i = p_value['nu_i']
    nu_g = p_value['nu_g']
    Em_v = p_value['Em_v']
    Em_i = p_value['Em_i']
    Em_g = p_value['Em_g']
    Eb_v_g = p_value['Eb_v_g']
    Eb_v_2g = p_value['Eb_v_2g']
    Eb_2g = p_value['Eb_2g']
    Ef_v = p_value['Ef_v']
    a0 = p_value['a0']
    Omega = p_value['Omega']
    f = p_value['f']
    b = p_value['b']
    k_B = p_value['k_B']
    gamma_b = p_value['gamma_b']
    B = p_value['B']
    r_ppt = p_value['r_ppt']
    d = p_value['d']
    N_ppt = p_value['N_ppt']
    rho = p_value['rho']
    Z_i = p_value['Z_i']

    # set print precision digits
    np.set_printoptions(precision=10)
    print('nu_v =', nu_v)
    print('nu_i =', nu_i)
    print('nu_g =', nu_g)
    print('Em_v =', Em_v)
    print('Em_i =', Em_i)
    print('Em_g =', Em_g)
    print('Eb_v_g =', Eb_v_g)
    print('Eb_v_2g =', Eb_v_2g)
    print('Eb_2g =', Eb_2g)
    print('Ef_v =', Ef_v)
    print('a0 =', a0)
    print('Omega =', Omega)
    print('f =', f)
    print('b =', b)
    print('k_B =', k_B)
    print('gamma_b =', gamma_b)
    print('B =', B)
    print('r_ppt =', r_ppt)
    print('d =', d)
    print('N_ppt =', N_ppt)
    print('rho =', rho)
    print('Z_i =', Z_i)
    
    return [str(nu_v), str(nu_i), str(nu_g), str(Em_v), str(Em_i), str(Em_g), str(Eb_v_g), str(Eb_v_2g), str(Eb_2g), str(Ef_v), str(a0), str(Omega), str(f), str(b), str(k_B), str(gamma_b), str(B), str(r_ppt), str(d), str(N_ppt), str(rho), str(Z_i)]

if __name__ == '__main__':
    
    list_of_parameters = read_data()
    
    # Run the program and capture the output
    result = subprocess.run(['./build/main'] + list_of_parameters, capture_output=True, text=True)
    
    # Check for errors
    if result.returncode != 0:
        print("Execution failed")
        print(result.stderr)
        exit(1)
    
    print(result.stdout)