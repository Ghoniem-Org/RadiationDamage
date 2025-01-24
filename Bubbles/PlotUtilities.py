import matplotlib.pyplot as plt
import os

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