import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

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

def parse_output(output):
    lines = output.strip().split('\n')
    data = [list(map(float, line.strip().split())) for line in lines]
    return np.array(data)

class FigureGenerator:
    def __init__(self, fig_width=8, fig_height=6, xlabel='X-axis', ylabel='Y-axis',
                 x_lim=None, y_lim=None, legend_fontsize=16, tick_fontsize=18,
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

        if os.path.exists(self.fig_path) is False:
            os.makedirs(self.fig_path)
            
        full_path = os.path.join(self.fig_path, save_filename)
        plt.savefig(full_path)
        if show:
            plt.show()
        plt.close()

def create_bubble_plots(t_eval, sol, t_span, fig_path, Omega, font1, font2):
    X_axis_lim = (10 * t_span[0], t_span[1])
    PointDefects = FigureGenerator(xlabel='Time [s]', ylabel='Concentration [$cm^{-3}$]',
                                   x_lim=X_axis_lim, y_lim=(1e6, 1e18), legend_fontsize=font1, tick_fontsize=font2, fig_path=fig_path, loglog=True)

    Gas = FigureGenerator(xlabel='Time [s]', ylabel='Atoms/ Bubble',
                          x_lim=X_axis_lim, legend_fontsize=font1, tick_fontsize=font2, fig_path=fig_path, logx=True)

    Radius = FigureGenerator(xlabel='Time [s]', ylabel='Avg. Radius [nm]',
                             x_lim=X_axis_lim, legend_fontsize=font1, tick_fontsize=font2, fig_path=fig_path, logx=True)  # Use logx=True for log-x, logy=True for log-y, or loglog=True for both

    data1 = [
        {'x': t_eval, 'y': sol[:, 0]/Omega/1e6, 'label': '$C_v$', 'color': 'blue'},
        {'x': t_eval, 'y': sol[:, 1]/Omega/1e6, 'label': '$C_i$', 'color': 'green'},
        {'x': t_eval, 'y': sol[:, 2]/Omega/1e6, 'label': '$C_g$', 'color': 'red'}
    ]

    data2 = [
        {'x': t_eval, 'y': sol[:, 3] / Omega / 1e6, 'label': '$C_{gv}$', 'color': 'blue'},
        {'x': t_eval, 'y': sol[:, 4] / Omega / 1e6, 'label': '$C_{g2v}$', 'color': 'green'},
        {'x': t_eval, 'y': sol[:, 5] / Omega / 1e6, 'label': '$C_{2g}$', 'color': 'red'}
    ]

    data3 = [
        {'x': t_eval, 'y': sol[:, 6] / Omega / 1e6, 'label': '$C^{*}$', 'color': 'blue'},
        {'x': t_eval, 'y': sol[:, 7] / Omega / 1e6, 'label': '$C_{B}$', 'color': 'red'},
    ]

    data4 = [
        {'x': t_eval, 'y': sol[:, 8], 'label': '$m$', 'color': 'blue'},
        {'x': t_eval, 'y': sol[:, 10], 'label': '$m_{ppt}$', 'color': 'red'},
    ]

    data5 = [
        {'x': t_eval, 'y': sol[:, 9]*1e9, 'label': '$R$', 'color': 'blue'},
        {'x': t_eval, 'y': sol[:, 11]*1e9, 'label': '$R_{ppt}$', 'color': 'red'},
    ]

    PointDefects.plot_data(data1, save_filename='point_defect.png')
    PointDefects.plot_data(data2, save_filename='clusters.png')
    PointDefects.plot_data(data3, save_filename='bubbles.png')
    Gas.plot_data(data4, save_filename='gas.png')
    Radius.plot_data(data5, save_filename='radius.png')
