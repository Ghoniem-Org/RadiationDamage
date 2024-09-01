import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import subprocess

from Utilities import read_data, parse_output
from PlotUtilities import create_bubble_plots

if __name__ == '__main__':
    
    # Read the data from the Excel file
    # list_of_parameters = read_data()
    parameters_dict = read_data()
    
    # Convert dictionary to list of keyword arguments
    keyword_arguments = []
    for key, value in parameters_dict.items():
        keyword_arguments.append(f'--{key}={value}')

    # Solve the ODE by running c++ executable
    result = subprocess.run(['./build/main'] + keyword_arguments, capture_output=True, text=True)
    
    # Print the output
    print(result.stdout)

    # Check for errors
    if result.returncode != 0:
        print("Execution failed")
        print(result.stderr)
        exit(1)
    
    sol = parse_output(result.stdout)
    print('Parsing complete ')

    # TODO: Perhalps pass time parameters as arguments to c++ executable?
    # Create bubble plots
    t_span = (1e-6, 1e6)
    # Time points at which to solve the ODE
    t_eval = np.logspace(np.log10(t_span[0]), np.log10(t_span[1]), 200)
    create_bubble_plots(t_eval, sol, t_span, fig_path = './figures', Omega=float(list_of_parameters[11]), font1=16, font2=18)
    print('Bubble plots created')