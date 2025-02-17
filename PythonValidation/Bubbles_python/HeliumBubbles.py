
import numpy as np
from scipy.integrate import odeint
from Derivatives import derivatives
from BubblePlots import create_plots
from Utilities import init_conditions
import os
print(os.getcwd())


y0, floor, a_length = init_conditions()
t_span = (1e-6, 1e6)

# Time points at which to solve the ODE
t_eval = np.logspace(np.log10(t_span[0]), np.log10(t_span[1]), 200)

# Solve the system of ODEs
sol = odeint(derivatives, y0, t_eval, atol=floor)
create_plots(t_eval, sol, t_span)
