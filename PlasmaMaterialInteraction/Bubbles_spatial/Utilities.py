import pandas as pd
import numpy as np

def read_data():
    # Load the Excel file
    df_steel = pd.read_excel('PlasmaMaterialInteraction\Bubbles_spatial\material_data.xlsx', sheet_name='steel', engine='openpyxl')

    # Convert the DataFrame to a dictionary
    # Assuming 'key' and 'value' are column names in your Excel file
    # parameters = df_steel.set_index('key').to_dict()['value']

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
    # np.set_printoptions(precision=10)
    # print('nu_v =', nu_v)
    # print('nu_i =', nu_i)
    # print('nu_g =', nu_g)
    # print('Em_v =', Em_v)
    # print('Em_i =', Em_i)
    # print('Em_g =', Em_g)
    # print('Eb_v_g =', Eb_v_g)
    # print('Eb_v_2g =', Eb_v_2g)
    # print('Eb_2g =', Eb_2g)
    # print('Ef_v =', Ef_v)
    # print('a0 =', a0)
    # print('Omega =', Omega)
    # print('f =', f)
    # print('b =', b)
    # print('k_B =', k_B)
    # print('gamma_b =', gamma_b)
    # print('B =', B)
    # print('r_ppt =', r_ppt)
    # print('d =', d)
    # print('N_ppt =', N_ppt)
    # print('rho =', rho)
    # print('Z_i =', Z_i)
    
    # return [str(nu_v), str(nu_i), str(nu_g), str(Em_v), str(Em_i), str(Em_g), str(Eb_v_g), str(Eb_v_2g), str(Eb_2g), str(Ef_v), str(a0), str(Omega), str(f), str(b), str(k_B), str(gamma_b), str(B), str(r_ppt), str(d), str(N_ppt), str(rho), str(Z_i)]
    return p_value

# def parse_output(output):
#     lines = output.strip().split('\n')
#     data = [list(map(float, line.strip().split())) for line in lines]
#     return np.array(data)

def parse_output_spatial(stdout_str, Nx, Ns, num_time_points, excel_filename="results_spatial.xlsx"):
    """
    Parse the stdout from the C++ code into a NumPy array and save results to an Excel file.

    Parameters
    ----------
    stdout_str : str
        The combined stdout from the C++ executable.
    Nx : int
        Number of spatial nodes.
    Ns : int
        Number of species.
    num_time_points : int
        Number of time points (including the initial condition).
    excel_filename : str
        Name of the Excel file to write sheets to.

    Returns
    -------
    sol : np.ndarray
        Array with shape (num_time_points, Ns, Nx).
    times : list of float
        List of time values corresponding to each slice in `sol`.
    """
    lines = stdout_str.strip().splitlines()
    idx_line = 0
    num_lines = len(lines)

    # Containers to accumulate time values and data blocks
    times = []
    data_blocks = []

    # Loop to parse output
    while idx_line < num_lines:
        line = lines[idx_line].strip()

        if line.startswith("Time:"):
            # Extract the timestamp
            parts = line.split()
            t_val = float(parts[1])  # Extract timestamp as float
            times.append(t_val)
            idx_line += 1  # Move to the next line after "Time:"

            # Read the next Nx lines as data for this timestamp
            block_lines = lines[idx_line : idx_line + Nx]
            idx_line += Nx

            block_matrix = []
            for bl in block_lines:
                cols = bl.strip().split()
                float_cols = list(map(float, cols))
                block_matrix.append(float_cols)
            block_matrix = np.array(block_matrix).T  # shape (Ns, Nx)
            data_blocks.append(block_matrix)

        else:
            # Skip lines not starting with "Time:" (e.g., header/debug info)
            idx_line += 1

    # Ensure we have the expected number of time points
    if len(data_blocks) != num_time_points:
        raise ValueError(
            f"Expected {num_time_points} time points, but parsed {len(data_blocks)}."
        )

    # Stack all data blocks into shape (num_time_points, Ns, Nx)
    sol = np.stack(data_blocks, axis=0)  # shape: (num_time_points, Ns, Nx)

    # Write results to Excel
    with pd.ExcelWriter(excel_filename) as writer:
        for i, t_val in enumerate(times):
            df = pd.DataFrame(data_blocks[i])
            sheet_name = f"{t_val:.6e}"  # Use the exact timestamp as the sheet name
            df.to_excel(writer, sheet_name=sheet_name, index=False, header=False)

    return sol, times

