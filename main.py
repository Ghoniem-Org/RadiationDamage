import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import subprocess

from Utilities import read_data, parse_output
from Utilities import create_creep_plots

if __name__ == '__main__':
    
    # Read the data from the Excel file
    list_of_parameters = read_data()
    
    # Solve the ODE by running c++ executable
    result = subprocess.run(['./build/main'] + list_of_parameters, capture_output=True, text=True)
    
    # Print the output
    print(result.stdout)

    # Check for errors
    if result.returncode != 0:
        print("Execution failed")
        print(result.stderr)
        exit(1)
    
    sol = parse_output(result.stdout)
    print('Parsing complete ')

    # TODO: Perhalps pass time parameters and yes_log as arguments to c++ executable?
    # Create bubble plots
    t_span = (1e-10, 4e7)
    # Time points at which to solve the ODE
    t_eval = np.linspace(1e-10, 4e7, 200)
    create_creep_plots(t_eval, sol, t_span, fig_path = './figures', font1=16, font2=18)
    print('Creep plots created')
