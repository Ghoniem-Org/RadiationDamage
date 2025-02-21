import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import subprocess

from Utilities import read_data, parse_output_spatial
# from PlotUtilities import create_bubble_plots

if __name__ == '__main__':

    # Read the data from the Excel file
    parameters_dict = read_data()
    
    # Convert dictionary to list of keyword arguments
    keyword_arguments = []
    for key, value in parameters_dict.items():
        keyword_arguments.append(f'--{key}={value}')

    # Read from parameters.txt and append to keyword_arguements
    with open("./parameters.txt", "r") as file:
        for line in file:
            key, value = line.strip().split(maxsplit=1)  # Split into name and value
            # Convert value to int or float
            value = float(value) if key != "time_points" and key != "spatial_nodes" else int(value)
            keyword_arguments.append(f'--{key}={value}')

            if key == "time_points":
                time_points = value
            if key == "spatial_nodes":
                spatial_nodes = value


    print(keyword_arguments)
    # Solve the ODE by running c++ executable
    result = subprocess.run(['./build/Release/main_spatial_from_curve'] + keyword_arguments, capture_output=True, text=True)
    # import subprocess
    # excutable_path = r'C:\Users\Owner\Documents\Repos\RadiationDamage\PlasmaMaterialInteraction\Bubbles_spatial\build\Debug\main_spatial_from_curve.exe'
    # result = subprocess.run([excutable_path] + keyword_arguments, capture_output=True, text=True)

    # Print output
    print(result.stdout)

    # Check for errors
    if result.returncode != 0:
        print("Execution failed")
        print(result.stderr)
        exit(1)
    
    sol, times = parse_output_spatial(result.stdout, spatial_nodes, 13, time_points, excel_filename="results_spatial.xlsx")
    print('Parsing complete ')

    # # Create bubble plots
    # t_eval = np.logspace(np.log10(t0), np.log10(tf), time_points)
    # create_bubble_plots(t_eval, sol, (t0,tf), fig_path = './figures', Omega=float(parameters_dict['Omega']), font1=16, font2=18)
    # print('Bubble plots created')