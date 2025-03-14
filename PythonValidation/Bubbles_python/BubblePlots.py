from Utilities import get_props, FigureGenerator

font1 = 16
font2 = 18
fig_path = r'.\Figures'

props = get_props()
Omega = props['Omega']


def create_plots(t_eval, sol, t_span):
    X_axis_lim = (10 * t_span[0], t_span[1])
    PointDefects = FigureGenerator(xlabel='Time [s]', ylabel='Concentration [$cm^{-3}$]',
                                   x_lim=X_axis_lim, y_lim=(1e6, 1e18), fig_path=fig_path, loglog=True)

    Gas = FigureGenerator(xlabel='Time [s]', ylabel='Atoms/ Bubble',
                          x_lim=X_axis_lim, fig_path=fig_path, logx=True)

    Radius = FigureGenerator(xlabel='Time [s]', ylabel='Avg. Radius [nm]',
                             x_lim=X_axis_lim, fig_path=fig_path, logx=True)  # Use logx=True for log-x, logy=True for log-y, or loglog=True for both

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
