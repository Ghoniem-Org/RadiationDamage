import matplotlib.pyplot as plt
import os

class FigureGenerator:
    def __init__(self, fig_width=8, fig_height=6, xlabel='X-axis', ylabel='Y-axis',
                 xticks=None, yticks=None,
                 x_lim=None, y_lim=None, legend_fontsize=16, tick_fontsize=18,
                 grid=False, legend_loc='best', fig_path='./', logx=False,
                 logy=False, loglog=False, linewidth=3, marker='o', facecolors='cyan'):
        self.fig_width = fig_width
        self.fig_height = fig_height
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xticks = xticks
        self.yticks = yticks
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
            
            if 'label' in dataset:

                plot_func(dataset['x'], dataset['y'], label=dataset['label'],
                        color=dataset.get('color', None),
                        linewidth=dataset.get('linewidth', 2))

                if 'scatter' in dataset:
                    plt.scatter(dataset['x'], dataset['y'], color=dataset.get('scatter_color', 'red'),
                                s=dataset.get('scatter_size', 50), label=dataset.get('scatter_label', None),
                                marker=dataset.get('marker', 'o'), facecolors=dataset.get('facecolors', 'none'))

            else:
                plot_func(dataset['x'], dataset['y'], color=dataset.get('color', None),
                        linewidth=dataset.get('linewidth', 2))
                
                if 'scatter' in dataset:
                    plt.scatter(dataset['x'], dataset['y'], color=dataset.get('scatter_color', 'red'),
                                s=dataset.get('scatter_size', 50),
                                marker=dataset.get('marker', 'o'), facecolors=dataset.get('facecolors', 'none'))
                    
        plt.xlabel(self.xlabel, fontsize=self.tick_fontsize)
        plt.ylabel(self.ylabel, fontsize=self.tick_fontsize)
        if self.x_lim:
            plt.xlim(self.x_lim)
        if self.y_lim:
            plt.ylim(self.y_lim)
        if self.xticks:
            plt.xticks(ticks=self.xticks, fontsize=self.tick_fontsize)
        else:
            plt.xticks(fontsize=self.tick_fontsize)
        if self.yticks:
            plt.yticks(ticks=self.yticks, fontsize=self.tick_fontsize)
        else:
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

def create_creep_plots(t_eval, sol, t_span, fig_path, font1, font2):
    # X_axis_lim = (10 * t_span[0], t_span[1])
    mobile_dislocation_density = FigureGenerator(xlabel='Time [s]', ylabel='Mobile Dislocation Density [1/m^2]',
                                                 yticks=[0e15, 0.5e15, 1e15, 1.5e15, 2e15, 2.5e15, 3e15],
                                   x_lim=(0, 4e7), y_lim=(0, 3e15), legend_fontsize=font1, tick_fontsize=font2, fig_path=fig_path)

    static_dislocation_density = FigureGenerator(xlabel='Time [s]', ylabel='Static Dislocation Density [1/m^2]',
                                                 yticks=[0e12 ,0.5e12, 1e12, 1.5e12, 2e12, 2.5e12],
                          x_lim=(0, 4e7), y_lim=(0, 2.5e12), legend_fontsize=font1, tick_fontsize=font2, fig_path=fig_path)

    boundary_dislocation_density = FigureGenerator(xlabel='Time [s]', ylabel='Boundry Dislocation Density [1/m^2]',
                                                   yticks=[1e11, 1.5e11, 2e11, 2.5e11, 3e11, 3.5e11, 4e11, 4.5e11],
                          x_lim=(0, 4e7), y_lim=(1e11, 4.5e11), legend_fontsize=font1, tick_fontsize=font2, fig_path=fig_path)

    subgrain_radius = FigureGenerator(xlabel='Time [s]', ylabel='Subgrain Radius [m]',
                                      yticks=[1e-5, 1.5e-5, 2e-5, 2.5e-5, 3e-5],
                          x_lim=(0, 4e7), y_lim=(1e-5, 3e-5), legend_fontsize=font1, tick_fontsize=font2, fig_path=fig_path)

    strain = FigureGenerator(xlabel='Time [hours]', ylabel='Boundry Dislocation Density [1/m^2]',
                                yticks=[0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4],
                          legend_fontsize=font1, tick_fontsize=font2, fig_path=fig_path)
    

    data1 = [
        {'x': t_eval, 'y': sol[:, 0], 'color': 'blue'}
    ]

    data2 = [
        {'x': t_eval, 'y': sol[:, 1], 'color': 'blue'}
    ]

    data3 = [
        {'x': t_eval, 'y': sol[:, 2], 'color': 'blue'}
    ]

    data4 = [
        {'x': t_eval, 'y': sol[:, 3], 'color': 'red'}
    ]

    data5 = [
        {'x': t_eval / 3600 , 'y': sol[:, 4]*1e2, 'color': 'red'}
    ]

    mobile_dislocation_density.plot_data(data1, save_filename='mobile_dislocation_density.png')
    static_dislocation_density.plot_data(data2, save_filename='static_dislocation_density.png')
    boundary_dislocation_density.plot_data(data3, save_filename='boundary_dislocation_density.png')
    subgrain_radius.plot_data(data4, save_filename='subgrain_radius.png')
    strain.plot_data(data5, save_filename='strain.png')