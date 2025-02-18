import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import subprocess

from Utilities import read_data, parse_output_spatial
# from PlotUtilities import create_bubble_plots

if __name__ == '__main__':
    
    # 
    # Read the data from the Excel file
    parameters_dict = read_data()
    
    # Convert dictionary to list of keyword arguments
    keyword_arguments = []
    for key, value in parameters_dict.items():
        keyword_arguments.append(f'--{key}={value}')

    # Append time parameters
    t0 = 1e-6
    tf = 1e6
    time_points = 200
    keyword_arguments.append(f'--t0={t0}')
    keyword_arguments.append(f'--tf={tf}')
    keyword_arguments.append(f'--time_points={time_points}')

    # Solve the ODE by running c++ executable
    # result = subprocess.run(['./build/main_spatial'] + keyword_arguments, capture_output=True, text=True)
    import subprocess
    excutable_path = r'C:\Users\Owner\Documents\Repos\RadiationDamage\PlasmaMaterialInteraction\Bubbles_spatial\build\Debug\main_spatial_from_curve.exe'
    result = subprocess.run([excutable_path] + keyword_arguments, capture_output=True, text=True)

    # Print output
    print(result.stdout)

    # Check for errors
    if result.returncode != 0:
        print("Execution failed")
        print(result.stderr)
        exit(1)
    
    sol, times = parse_output_spatial(result.stdout, 20, 13, time_points, excel_filename="results_spatial.xlsx")
    print('Parsing complete ')

    # # Create bubble plots
    # t_eval = np.logspace(np.log10(t0), np.log10(tf), time_points)
    # create_bubble_plots(t_eval, sol, (t0,tf), fig_path = './figures', Omega=float(parameters_dict['Omega']), font1=16, font2=18)
    # print('Bubble plots created')