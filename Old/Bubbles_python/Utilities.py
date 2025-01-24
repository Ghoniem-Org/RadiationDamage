import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

font1 = 16
font2 = 18

# Input parameters
T = 625 + 273
G = 3e-3
he_2_dpa = 5e-6

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


def get_props():

    # Reaction frequencies
    alpha = 48 * nu_i * np.exp(-Em_i / (k_B * T))
    beta = 48 * nu_g * np.exp(-Em_g / (k_B * T))
    gamma = 48 * nu_v * np.exp(-Em_v / (k_B * T))

    # Thermal emission frequencies
    e1 = np.exp(-Eb_v_g / (k_B * T))
    e2 = np.exp(-Eb_v_2g / (k_B * T))
    e4 = np.exp(-Ef_v / (k_B * T))
    e5 = np.exp(-Eb_2g / (k_B * T))

    # Diffusion coefficients
    D_i = (a0 ** 2 / 48) * alpha
    D_v = (a0 ** 2 / 48) * gamma
    D_g = (a0 ** 2 / 48) * beta

    # Equilibrium vacancy concentration
    C_v_e = np.exp(-Ef_v / (k_B * T))
    C_ppt = N_ppt * Omega

    # Production terms
    P = f * G
    delta = b * G
    G_He = he_2_dpa * G

    return {'alpha': alpha, 'beta': beta, 'gamma': gamma, 'e1': e1,
            'e2': e2, 'e4': e4, 'e5': e5, 'delta': delta, 'D_i': D_i,
            'D_v': D_v, 'D_g': D_g, 'P': P, 'G_He': G_He, 'C_v_e': C_v_e,
            'Omega': Omega, 'C_ppt': C_ppt, 'rho': rho, 'a0': a0, 'd': d,
            'Z_i': Z_i, 'r_ppt': r_ppt}


def get_bubble_props(R, m):
    pressure = m * k_B * T / (4. * np.pi * R ** 3 / 3 - m * B)
    work = (2 * gamma_b / R - pressure) * Omega
    Eb_v_B = Ef_v + work
    e3 = np.exp(-Eb_v_B / (k_B * T))
    epsilon = (4 * np.pi/48) * (R / a0)
    return e3, epsilon, work


def init_conditions():
    floor = 1e-20
    array_length = 13
    y0 = np.full(array_length, floor)
    y0[8] = 2
    y0[9] = 5e-10
    y0[10] = 2
    y0[11] = 5e-10
    return y0, floor, array_length


class FigureGenerator:
    def __init__(self, fig_width=8, fig_height=6, xlabel='X-axis', ylabel='Y-axis',
                 x_lim=None, y_lim=None, legend_fontsize=font1, tick_fontsize=font2,
                 grid=False, legend_loc='best', fig_path='./', logx=False,
                 logy=False, loglog=False, linewidth=3, marker='o', facecolors='cyan'):
        self.fig_width = fig_width
        self.fig_height = fig_height
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.x_lim = x_lim
        self.y_lim = y_lim
        self.legend_fontsize = legend_fontsize
        self.tick_fontsize = tick_fontsize
        self.grid = grid
        self.legend_loc = legend_loc
        self.fig_path = fig_path
        self.logx = logx
        self.logy = logy
        self.loglog = loglog
        self.linewidth = linewidth

    def plot_data(self, data, save_filename='figure.png', show=False):
        plt.figure(figsize=(self.fig_width, self.fig_height))

        # Determine the type of plot based on log scale options
        plot_func = plt.plot
        if self.loglog:
            plot_func = plt.loglog
        elif self.logx:
            plt.xscale('log')
        elif self.logy:
            plt.yscale('log')

        for dataset in data:
            plot_func(dataset['x'], dataset['y'], label=dataset['label'],
                      color=dataset.get('color', None),
                      linewidth=dataset.get('linewidth', 2))

            if 'scatter' in dataset:
                plt.scatter(dataset['x'], dataset['y'], color=dataset.get('scatter_color', 'red'),
                            s=dataset.get('scatter_size', 50), label=dataset.get('scatter_label', None),
                            marker=dataset.get('marker', 'o'), facecolors=dataset.get('facecolors', 'none'))

        plt.xlabel(self.xlabel, fontsize=self.tick_fontsize)
        plt.ylabel(self.ylabel, fontsize=self.tick_fontsize)
        if self.x_lim:
            plt.xlim(self.x_lim)
        if self.y_lim:
            plt.ylim(self.y_lim)
        plt.xticks(fontsize=self.tick_fontsize)
        plt.yticks(fontsize=self.tick_fontsize)
        plt.grid(self.grid)
        plt.legend(loc=self.legend_loc, fontsize=self.legend_fontsize)

        full_path = os.path.join(self.fig_path, save_filename)
        plt.savefig(full_path)
        if show:
            plt.show()
        plt.close()

