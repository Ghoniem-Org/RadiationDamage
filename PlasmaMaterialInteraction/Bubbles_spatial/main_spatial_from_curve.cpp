// #include "parameters_spatial.h"
#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <cvodes/cvodes.h>               // CVODE functions
#include <nvector/nvector_serial.h>      // N_Vector
#include <sundials/sundials_dense.h>     // dense solver
#include <sunlinsol/sunlinsol_dense.h>   // dense linear solver
#include <sunmatrix/sunmatrix_dense.h>   // dense SUNMatrix
#include <chrono>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Inline indexing function
inline int idx(int i, int j, int Ns, int Nx) {
    return j*Ns + i;
}

// // Declare and define file-level constants
// const double T = 625 + 273;
// // const double T = 300;
// const double G = 3e-3;
// const double he_2_dpa = 5e-6;

// Forward declarations for functions to get properties
std::map<std::string, double> get_properties(int argc, char *argv[]) {
// std::map<std::string, double> get_properties() {

    std::map<std::string, double> props;

    // // Loop through all arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.substr(0, 2) == "--") { // Check if argument starts with "--"
            auto delimiter_pos = arg.find('=');
            if (delimiter_pos != std::string::npos) {
                std::string key = arg.substr(2, delimiter_pos - 2); // Extract the key (without "--")
                std::string value_str = arg.substr(delimiter_pos + 1); // Extract the value as a string

                // Convert string to double
                try {
                    double value = std::stod(value_str);
                    props[key] = value; // Store in the map
                } catch (const std::invalid_argument& e) {
                    std::cerr << "Invalid argument for " << key << ": " << value_str << std::endl;
                    exit(1);
                } catch (const std::out_of_range& e) {
                    std::cerr << "Value out of range for " << key << ": " << value_str << std::endl;
                    exit(1);
                }
            } else {
                std::cerr << "Invalid argument format: " << arg << std::endl;
                exit(1);
            }
        } else {
            std::cerr << "Unexpected argument: " << arg << std::endl;
            exit(1);
        }
    }

    // Comment when using command line arguments
    // props["nu_v"] = 5000000000000;
    // props["nu_i"] = 50000000000000;
    // props["nu_g"] = 50000000000000;
    // props["Em_v"] = 1.4;
    // props["Em_i"] = 0.2;
    // props["Em_g"] = 0.2;
    // props["Eb_v_g"] = 2.4;
    // props["Eb_v_2g"] = 3.5;
    // props["Eb_2g"] = 0.79;
    // props["Ef_v"] = 1.6;
    // props["a0"] = 3.63e-10;
    // props["Omega"] = 1.195803675E-029;
    // props["f"] = 1e-1;
    // props["b"] = 10;
    // props["k_B"] = 8.617e-5;
    // props["gamma_b"] = 6.24e18;
    // props["B"] = 1.75e-29;
    // props["r_ppt"] = 1e-8;
    // props["d"] = 3e-5;
    // props["N_ppt"] = 1e16;
    // props["rho"] = 300000000000000;
    // props["Z_i"] = 1.2;
    
    // props["t0"] = 1e-6;
    // props["tf"] = 1e6;
    // props["time_points"] = 200;

    // Retrieve parameters from props map or use default values
    double nu_v = props["nu_v"];
    double nu_i = props["nu_i"];
    double nu_g = props["nu_g"];
    double Em_v = props["Em_v"];
    double Em_i = props["Em_i"];
    double Em_g = props["Em_g"];
    double Eb_v_g = props["Eb_v_g"];
    double Eb_v_2g = props["Eb_v_2g"];
    double Eb_2g = props["Eb_2g"];
    double Ef_v = props["Ef_v"];
    double a0 = props["a0"];
    double Omega = props["Omega"];
    double f = props["f"];
    double b = props["b"];
    double k_B = props["k_B"];
    double gamma_b = props["gamma_b"];
    double B = props["B"];
    double r_ppt = props["r_ppt"];
    double d = props["d"];
    double N_ppt = props["N_ppt"];
    double rho = props["rho"];
    double Z_i = props["Z_i"];

    // Parameters from parameters.txt
    double t0 = props["t0"];
    double tf = props["tf"];
    double time_points = props["time_points"];
    double spatial_nodes = props["spatial_nodes"];
    double xN = props["xN"];
    double A_P = props["A_P"];
    double B_P = props["B_P"];
    double C_P = props["C_P"];
    double A_G = props["A_G"];
    double B_G = props["B_G"];
    double C_G = props["C_G"];
    double T = props["T"];
    double G = props["G"];
    double he_2_dpa = props["he_2_dpa"];

    // Reaction frequencies
    double alpha = 48 * nu_i * exp(-Em_i / (k_B * T));
    double beta = 48 * nu_g * exp(-Em_g / (k_B * T));
    double gamma = 48 * nu_v * exp(-Em_v / (k_B * T));

    // double alpha = 48 * nu_i * exp(-Em_i / k_B * T);
    // double beta = 48 * nu_g * exp(-Em_g / k_B * T);
    // double gamma = 48 * nu_v * exp(-Em_v / k_B * T);
    // Thermal emission frequencies
    double e1 = exp(-Eb_v_g / (k_B * T));
    double e2 = exp(-Eb_v_2g / (k_B * T));
    double e4 = exp(-Ef_v / (k_B * T));
    double e5 = exp(-Eb_2g / (k_B * T));

    // Diffusion coefficients
    double D_i = (a0 * a0 / 48) * alpha;
    double D_v = (a0 * a0 / 48) * gamma;
    double D_g = (a0 * a0 / 48) * beta;

    // Equilibrium vacancy concentration
    double C_v_e = exp(-Ef_v / (k_B * T));
    double C_ppt = N_ppt * Omega;

    // Production terms
    // double P = f * G; //Fuction of space
    double delta = b * G;
    // double G_He = he_2_dpa * G; //Fuction of space

    // Populate the map with calculated properties
    props["alpha"] = alpha;
    props["beta"] = beta;
    props["gamma"] = gamma;
    props["e1"] = e1;
    props["e2"] = e2;
    props["e4"] = e4;
    props["e5"] = e5;
    props["delta"] = delta;
    props["D_i"] = D_i;
    props["D_v"] = D_v;
    props["D_g"] = D_g;
    // props["P"] = P;
    // props["G_He"] = G_He;
    props["C_v_e"] = C_v_e;
    props["C_ppt"] = C_ppt;
    props["rho"] = rho;
    props["a0"] = a0;
    props["d"] = d;
    props["Z_i"] = Z_i;
    props["r_ppt"] = r_ppt;
    props["k_B"] = k_B;
    props["gamma_b"] = gamma_b;
    props["B"] = B;
    props["T"] = T;
    props["Ef_v"] = Ef_v;
    props["t0"] = t0;
    props["tf"] = tf;
    props["time_points"] = time_points;
    props["spatial_nodes"] = spatial_nodes;
    props["xN"] = xN;
    props["A_P"] = A_P;
    props["B_P"] = B_P;
    props["C_P"] = C_P;
    props["A_G"] = A_G;
    props["B_G"] = B_G;
    props["C_G"] = C_G;

    return props;
}


std::tuple<double,double,double> get_bubble_props(double R, double m, std::map<std::string,double>& props){
    double k_B=props["k_B"];
    double T=props["T"];
    double B=props["B"];
    double gamma_b=props["gamma_b"];
    double Omega=props["Omega"];
    double Ef_v=props["Ef_v"];

    double pi=M_PI;
    double pressure = m*k_B*T/((4.0/3.0)*pi*R*R*R - m*B);
    double work=(2.0*gamma_b/R - pressure)*Omega;
    double Eb_v_B=Ef_v+work;
    double e3=exp(-Eb_v_B/(k_B*T));
    double epsilon=(4.0*pi/48.0)*(R/props["a0"]);

    return std::make_tuple(e3, epsilon, work);
}

// Compute production and loss terms, can be modified as needed
// For demonstration, we just reuse the logic from the original code, 
// now applied at each spatial node. If spatial dependence is required, 
// you could incorporate `x` into these functions.
static double compute_P(/* species, x, etc. */) {
    // Placeholder: return some P value if needed
    // In the original code, P was f*G (constant).
    // If you want position-dependent P, incorporate x.
    return 0.0; 
}

static double compute_L(/* species, x, etc. */) {
    // Placeholder: return some L value if needed
    return 0.0; 
}

// User data structure
// changed here:
struct Parameters {
    std::map<std::string, double> props;
    int Nx;
    int Ns;
    double x0;
    double xN;
    std::vector<double> P;
    std::vector<double> G_He;
};


// RHS function
int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data) {
    Parameters* params = static_cast<Parameters*>(user_data);
    std::map<std::string,double>& props = params->props;

    int Nx = params->Nx; 
    int Ns = params->Ns;
    double x0 = params->x0;
    double xN = params->xN;
    double dx = (xN - x0)/(Nx - 1);

    // Extract parameters (same as original code)
    double alpha = props["alpha"];
    double beta = props["beta"];
    double gamma_val = props["gamma"];
    double e1 = props["e1"];
    double e2 = props["e2"];
    double e4 = props["e4"];
    double e5 = props["e5"];
    double delta = props["delta"];
    // Changed here
    const std::vector<double>& P = params->P;
    const std::vector<double>& G_He = params->G_He;
    double Omega = props["Omega"];
    double C_ppt = props["C_ppt"];
    double rho = props["rho"];
    double a0 = props["a0"];
    double d = props["d"];
    double Z_i = props["Z_i"];
    double r_ppt = props["r_ppt"];
    double Ef_v = props["Ef_v"];
    double gamma_b = props["gamma_b"];
    double B = props["B"];
    double k_B = props["k_B"];
    double T = props["T"];

    // We now have 13 species at each point: 
    // Let's assume the order of species is the same as original:
    // 0: C_v
    // 1: C_i
    // 2: C_g
    // 3: C_gv
    // 4: C_2gv
    // 5: C_2g
    // 6: C_star
    // 7: C_b
    // 8: m
    // 9: R
    // 10: m_ppt
    // 11: R_pt
    // 12: M_gb

    // If you have species-specific diffusion coefficients, define them:
    // For demonstration, assume all species have the same diffusion coefficient set 
    // or retrieve from props if available. The original code had D_i, D_v, D_g...
    // Let's assume C_v uses D_v, C_i uses D_i, C_g uses D_g, and others have negligible diffusion (or set to 0).
    // Modify as needed.
    double D_v = props["D_v"]; // for vacancies
    double D_i = props["D_i"]; // for interstitials
    double D_g = props["D_g"]; // for helium maybe
    // For other species, set D=0 if they don't diffuse or if not defined.
    // This is problem-dependent. We will just demonstrate a scenario:
    // double D_array[13];
    // Assign diffusion coefficients to species:
    // This is an assumption:
    // D_array[0] = D_v; // C_v
    // D_array[1] = D_i; // C_i
    // D_array[2] = D_g; // C_g
    // // For others, assume no diffusion:
    // for (int i = 3; i < 13; i++) D_array[i] = 0.0;

    // Prevents time evolution at boundaries
    for (int i = 0; i < Ns; i++) {
        // NV_Ith_S(y, idx(i,0,Ns,Nx)) = 0.0;
        // NV_Ith_S(y, idx(i,Nx - 1,Ns,Nx)) = 0.0;
        NV_Ith_S(ydot, idx(i,0,Ns,Nx)) = 0.0;
        NV_Ith_S(ydot, idx(i,Nx - 1,Ns,Nx)) = 0.0;
    }

    // Compute derivatives
    for (int j = 1; j < Nx - 1; j++) {
        double x = x0 + j*dx;

        // Retrieve fields at node j:
        // For convenience, de fine some shorthand variables per node:
        // We'll directly use NV_Ith_S(y, idx(...)) in equations as before, 
        // but now for each j we have them spatially resolved.

        // We'll follow the same ODE structure as the original single-point code, 
        // but apply it at each node. The original code computed everything for a single cell;
        // now we do it for each j.

        // Retrieve needed species at node j:
        double C_v_val   = NV_Ith_S(y, idx(0,j,Ns,Nx));
        double C_i_val   = NV_Ith_S(y, idx(1,j,Ns,Nx));
        double C_g_val   = NV_Ith_S(y, idx(2,j,Ns,Nx));
        double C_gv_val  = NV_Ith_S(y, idx(3,j,Ns,Nx));
        double C_2gv_val = NV_Ith_S(y, idx(4,j,Ns,Nx));
        double C_2g_val  = NV_Ith_S(y, idx(5,j,Ns,Nx));
        double C_star_val= NV_Ith_S(y, idx(6,j,Ns,Nx));
        double C_b_val   = NV_Ith_S(y, idx(7,j,Ns,Nx));
        double m_val     = NV_Ith_S(y, idx(8,j,Ns,Nx));
        double R_val     = NV_Ith_S(y, idx(9,j,Ns,Nx));
        double m_ppt_val = NV_Ith_S(y, idx(10,j,Ns,Nx));
        double R_pt_val  = NV_Ith_S(y, idx(11,j,Ns,Nx));
        double M_gb_val  = NV_Ith_S(y, idx(12,j,Ns,Nx));

        // Bubble properties at this node:
        double e3, epsilon, work;
        std::tie(e3, epsilon, work) = get_bubble_props(R_val, m_val, props);
        double e3_prime, epsilon_ppt, work_ppt;
        std::tie(e3_prime, epsilon_ppt, work_ppt) = get_bubble_props(R_pt_val, m_ppt_val, props);

        double k_square = 4 * M_PI * R_val * C_b_val / Omega + rho;
        double C_gb_val_local = a0 * a0 * std::sqrt(k_square)/(8*d);

        double bubble_sink = 4 * M_PI * R_val * C_b_val / Omega;
        double precipitate_sink = 4 * M_PI * R_pt_val * C_ppt / Omega;

        double C_s_v = (a0*a0/48.0)*(rho + bubble_sink);
        double C_s_i = (a0*a0/48.0)*(Z_i * rho + bubble_sink);

        // For diffusion, we need neighbors:
        // Handle boundaries by no-flux mirroring:
        for (int i = 0; i < Ns; i++) {
            double C_left, C_center, C_right;
            // if (j == 0 || j == Nx-1){
            //     NV_Ith_S(y, idx(i,j,Ns,Nx)) = 0.0
            // }
            
            // if (j != 0 || j != Nx-1){
            C_left   = NV_Ith_S(y, idx(i,j-1,Ns,Nx));
            C_center = NV_Ith_S(y, idx(i,j,Ns,Nx));
            C_right  = NV_Ith_S(y, idx(i,j+1,Ns,Nx));
            double d2Cdx2 = (C_right - 2*C_center + C_left)/(dx*dx);
            
            double val = 0.0; // placeholder for ydot

            if (i == 0) {
                // dC_v/dt
                val = P[j] + (beta*e1+delta)*C_gv_val 
                      - (alpha*C_i_val + beta*C_g_val + gamma_val*(C_s_v + C_gv_val + 2*(C_2g_val + C_2gv_val) + 3*C_star_val))*C_v_val;
                // Add diffusion:
                val += D_v*d2Cdx2;
            } else if (i == 1) {
                // dC_i/dt
                val = P[j] - alpha*(C_v_val + C_gv_val + 2*C_2gv_val + 3*C_star_val + C_s_i)*C_i_val;
                val += D_i*d2Cdx2;
            } else if (i == 2) {
                // dC_g/dt
                val = (G_He[j] - beta*C_g_val*(C_v_val + 2*C_g_val + C_gv_val + 2*C_2g_val + 2*C_2gv_val + C_gb_val_local + epsilon*C_b_val)
                       + delta*(C_gv_val + 2*C_2gv_val + 4*C_2g_val + 3*C_star_val + m_val*C_b_val + M_gb_val + m_ppt_val*C_ppt)
                       + alpha*C_i_val*C_gv_val + beta*(e1*C_gv_val + e2*C_2gv_val));
                val += D_g*d2Cdx2;
            } else if (i == 3) {
                // dC_gv/dt
                val = beta*C_g_val*C_v_val + (beta*e2+2*delta)*C_2gv_val
                      - (beta*e1 + delta + alpha*C_i_val + beta*C_g_val)*C_gv_val;
            } else if (i == 4) {
                // dC_2gv/dt
                val = beta*C_g_val*C_gv_val + 3*delta*C_star_val + 2*gamma_val*C_v_val*C_2g_val
                      - (2*beta*C_g_val + 2*delta + beta*e2 + 2*alpha*C_i_val)*C_2gv_val;
            } else if (i == 5) {
                // dC_2g/dt
                val = alpha*C_i_val*C_2gv_val + 2*beta*C_g_val*C_g_val
                      - (2*delta + 2*gamma_val*C_v_val + 2*beta*C_g_val)*C_2g_val;
            } else if (i == 6) {
                // dC_star/dt
                val = 2*beta*C_g_val*(C_2gv_val + C_2g_val)
                      - 3*C_star_val*(delta + alpha*C_i_val + beta*C_g_val + gamma_val*C_v_val);
            } else if (i == 7) {
                // dC_b/dt
                val = ((12*beta*C_g_val + 9*gamma_val*C_v_val)*C_star_val)/m_val;
            } else if (i == 8) {
                // dm/dt
                val = epsilon*beta*C_g_val - delta*m_val;
            } else if (i == 9) {
                // dR/dt
                // If negative, set to zero
                {
                    double tmp = (a0*a0/R_val)*(gamma_val*C_v_val - alpha*C_i_val - gamma_val*(e3 - e4));
                    if (tmp < 0) tmp = 0.0;
                    val = tmp;
                }
            } else if (i == 10) {
                // dm_ppt/dt
                double epsilon_ppt_local = epsilon_ppt; // from bubble props
                val = epsilon_ppt_local*beta*C_g_val - delta*m_ppt_val;
            } else if (i == 11) {
                // dR_pt/dt
                {
                    double r_p_equiv = std::sqrt(R_pt_val*R_pt_val + r_ppt*r_ppt);
                    double tmp = (a0*a0 / r_p_equiv)*(gamma_val*C_v_val - alpha*C_i_val - gamma_val*(e3_prime - e4));
                    if (tmp < 0) tmp = 0.0;
                    val = tmp;
                }
            } else if (i == 12) {
                // dM_gb/dt
                val = beta*C_gb_val_local*C_g_val - delta*M_gb_val;
            }

            NV_Ith_S(ydot, idx(i,j,Ns,Nx)) = val;
        } // end species loop
    } // end spatial loop

    return 0;
}

// Main function similar to original, but with Nx, Ns adjustments
int main(int argc, char *argv[]) {
    auto props = get_properties(argc, argv);
    // auto props = get_properties();

    Parameters params;
    params.props = props;
    params.Ns = 13; // 13 species
    params.Nx = int(props["spatial_nodes"]); // 100 spatial nodes
    params.x0 = 0.0;
    params.xN = props["xN"];

    params.P.resize(params.Nx);
    params.G_He.resize(params.Nx);

    // Changed here
    const double x_max = params.xN;  // 1000 Å in meters
    const double dx = x_max / (params.Nx - 1);  // Spatial step size

    const double A_P = props["A_P"];
    const double B_P = props["B_P"];
    const double C_P = props["C_P"];

    const double A_G = props["A_G"];
    const double B_G = props["B_G"];
    const double C_G = props["C_G"];

    for (int j = 0; j < params.Nx; j++) {
        double x = j * dx;  // Convert index to physical distance
        params.P[j] = A_P * pow(x, B_P) * exp(-C_P * x);
        params.G_He[j] = A_G * pow(x, B_G) * exp(-C_G * x);
    }

    // ODE system size
    sunindextype Nx = params.Nx;
    sunindextype Ns = params.Ns;
    sunindextype neq = Nx * Ns;

    double t0 = props["t0"];
    double tf = props["tf"];
    double reltol = 1e-8;
    double abstol = 1e-20;

   

    SUNContext sunctx;
    SUNContext_Create(NULL, &sunctx);

    N_Vector y = N_VNew_Serial(neq, sunctx);
    for (int j = 0; j < Nx; j++) {
        double x = params.x0 + j*(params.xN - params.x0)/(Nx-1);
        for (int i = 0; i < Ns; i++) {
            // Set initial conditions
            // Just as in original code, you might set them from props or some initial guess
            // For now, set a small floor:
            if (i == 8 || i == 10) NV_Ith_S(y, idx(i,j,Ns,Nx)) = 2.0;
            else if (i == 9 || i == 11) NV_Ith_S(y, idx(i,j,Ns,Nx)) = 5e-10;
            else NV_Ith_S(y, idx(i,j,Ns,Nx)) = 1e-20;
        }
    }

    //Dirichlet BC
    for (int i = 0; i < Ns - 5; i++) {
        NV_Ith_S(y, idx(i, 0, Ns, Nx)) = 0;      // Left boundary (x = 0)
        NV_Ith_S(y, idx(i, Nx-1, Ns, Nx)) = 0;  // Right boundary (x = L)
    }

    void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
    CVodeSetUserData(cvode_mem, &params);
    CVodeInit(cvode_mem, rhs, t0, y);
    CVodeSetMaxNumSteps(cvode_mem, 100000);
    CVodeSetStabLimDet(cvode_mem, SUNTRUE);
    CVodeSetMinStep(cvode_mem, 1e-20);
    CVodeSStolerances(cvode_mem, reltol, abstol);

    SUNMatrix A = SUNDenseMatrix(neq, neq, sunctx);
    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
    CVodeSetLinearSolver(cvode_mem, LS, A);

    int time_points = (int)props["time_points"];
    std::vector<double> t_eval(time_points);
    double log_t_span_start = std::log10(t0);
    double log_t_span_end = std::log10(tf);
    double step = (log_t_span_end - log_t_span_start)/(time_points-1);
    for (int i = 0; i < time_points; i++) {
        t_eval[i] = std::pow(10, log_t_span_start + i*step);
    }

    double t = t0;
    
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "Time: " << t << std::endl;

    // Create a number of species x spatial nodes matrix
    std::vector<std::vector<double>> state_matrix(Nx, std::vector<double>(Ns));
    
    // Push initial
    for (int j = 0; j < Nx; j++) {
        for (int i = 0; i < Ns; i++) {
            state_matrix[j][i] = NV_Ith_S(y, idx(i,j,Ns,Nx));
            std::cout << state_matrix[j][i] << " ";
        }
        std::cout << "\n";
    }
    
    for (int i = 1; i < time_points; i++) {
        double tout = t_eval[i];
        int flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        if (flag < 0) {
            std::cerr << "Error in CVode at time " << t << " with flag " << flag << std::endl;
            break;
        }
        // To remove: print time
        std::cout << "Time: " << t << std::endl;

        // Create a number of species x spatial nodes matrix
        std::vector<std::vector<double>> state_matrix(Nx, std::vector<double>(Ns));
        
        for (int j = 0; j < Nx; j++) {
            for (int i = 0; i < Ns; i++) {
                state_matrix[j][i] = NV_Ith_S(y, idx(i,j,Ns,Nx));
                std::cout << state_matrix[j][i] << " ";
            }
            std::cout << "\n";
        }
        // results.push_back(state_matrix);
    }

    std::cout << std::scientific << std::setprecision(10);
    // for (auto &row : results) {
    //     for (auto &val : row) {
    //         std::cout << val << " ";
    //     }
    //     std::cout << "\n";
    // }

    N_VDestroy(y);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    CVodeFree(&cvode_mem);
    SUNContext_Free(&sunctx);
    return 0;
}
