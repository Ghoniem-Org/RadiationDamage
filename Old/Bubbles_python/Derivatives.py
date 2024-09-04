from Utilities import get_props, get_bubble_props, init_conditions
import numpy as np

# # # Read in basic parameters
# (T, nu_v, nu_g, nu_i, a0, b, G, appm_He, f, Em_g, Em_v, Em_i, Eb_v_g,
#  Eb_v_2g, Eb_2g, Ef_v, Z_i, rho, k_B, k_B_SI, gamma_b, B, r_ppt, d, N_ppt) = params.values()

# Calculate derived properties
props = get_props()
P = props['P']
alpha = props['alpha']
beta = props['beta']
gamma = props['gamma']
delta = props['delta']
e1 = props['e1']
e2 = props['e2']
e4 = props['e4']
e5 = props['e5']
Omega = props['Omega']
G_He = props['G_He']
C_ppt = props['C_ppt']
rho = props['rho']
a0 = props['a0']
d = props['d']
Z_i = props['Z_i']
r_ppt = props['r_ppt']

# Initialize the derivative array
y0, floor, a_length = init_conditions()

dydt = np.zeros(a_length)


# Define the function that returns the time derivatives of the ODEs
def derivatives(y, t):
    #
    # if y[9] <= y0[9]:
    #     y[9] = y0[9]
    # elif y[11] <= y0[11]:
    #     y[11] = y0[11]

    y = np.array(y)  # ensure y is a numpy array

  #  y_max = np.maximum(floor, y)

    (C_v, C_i, C_g, C_gv, C_2gv, C_2g, C_star, C_b, m, R, m_ppt, R_pt, M_gb) = y

    # Prepare time-dependent bubble-related parameters:
    e3, epsilon, work = get_bubble_props(R, m)
    e3_prime, epsilon_ppt, work_ppt = get_bubble_props(R_pt, m_ppt)
    k_square = 4*np.pi*R*C_b/Omega + rho
    C_gb = a0**2*np.sqrt(k_square)/(8*d)
    # Sink concentrations
    bubble_sink = 4*np.pi*R*C_b/Omega
    precipitate_sink = 4*np.pi*R_pt*C_ppt/Omega

    # counter = 0
    # if t < 1e-3 and counter <= 1:
    #     print('bubble_sink= ', format(bubble_sink, '.3e'))
    #     print('precipitate_sink= ', format(precipitate_sink, '.3e'))
    #     print('work=', format(work, '.3e'))
    #     counter += 1

    C_s_v = (a0 ** 2 / 48) * (rho + bubble_sink)
    C_s_i = (a0 ** 2 / 48) * (Z_i * rho + bubble_sink)

    #  Write down rate equations:
    #
    #  (dC_v/dt)
    dydt[0] = (P + (beta * e1 + delta) * C_gv
               - (alpha * C_i + beta * C_g
                  + gamma * (C_s_v + C_gv + 2 * (C_2g + C_2gv) + 3 * C_star)) * C_v)

    #  (dC_i/dt)
    dydt[1] = P - alpha * (C_v + C_gv + 2 * C_2gv + 3 * C_star + C_s_i) * C_i

    #  (dC_g/dt)
    dydt[2] = ((G_He - beta * C_g * (C_v + 2 * C_g + C_gv + 2 * C_2g + 2 * C_2gv
                                     + C_gb + epsilon * C_b)
               + delta * (C_gv + 2 * C_2gv + 4 * C_2g + 3 * C_star + m * C_b
                          + M_gb + m_ppt * C_ppt))
               + alpha * C_i * (C_gv) + beta * (e1 * C_gv + e2 * C_2gv))

    #  (dC_gv/dt)
    dydt[3] = (beta * C_g * C_v + (beta * e2 + 2 * delta) * C_2gv
               - (beta * e1 + delta + alpha * C_i + beta * C_g) * C_gv)

    #  (dC_2gv/dt)
    dydt[4] = (beta * C_g * C_gv + 3 * delta * C_star + 2 * gamma * C_v * C_2g
               - (2 * beta * C_g + 2 * delta + beta * e2
                  + 2 * alpha * C_i) * C_2gv)

    #  (dC_2g/dt)
    dydt[5] = (alpha * C_i * C_2gv + 2 * beta * C_g ** 2
               - (2 * delta + 2 * gamma * C_v + 2 * beta * C_g) * C_2g)

    #  (dC_star/dt)
    dydt[6] = (2 * beta * C_g * (C_2gv + C_2g)
               - 3 * C_star * (delta + alpha * C_i
                               + beta * C_g + gamma * C_v))

    #  (dC_b/dt)
    dydt[7] = ((12 * beta * C_g + 9 * gamma * C_v) * C_star) / m

    #  (dm/dt)
    dydt[8] = epsilon * beta * C_g - delta * m

    #  (dR/dt)
    dydt[9] = (a0**2/R) * (gamma * C_v - alpha * C_i - gamma * (e3 - e4))
    if dydt[9] < 0:
        dydt[9] = 0

    #  (dm_ppt/dt)
    dydt[10] = epsilon_ppt * beta * C_g - delta * m_ppt

    #  (dR_pt/dt)
    r_p_equiv = np.sqrt(R_pt**2+r_ppt**2)
    dydt[11] = (a0 ** 2 / r_p_equiv) * (gamma * C_v - alpha * C_i - gamma * (e3_prime - e4))
    if dydt[11] < 0:
        dydt[11] = 0

    #  (dM_gb/dt)
    dydt[12] = beta * C_gb * C_g - delta * M_gb

    return dydt
