#include "parameters.h"
#include <iomanip>
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <cvodes/cvodes.h>               // prototypes for CVODE functions
#include <nvector/nvector_serial.h>      // serial N_Vector types, functions, and macros
#include <sundials/sundials_dense.h>     // generic DENSE solver
#include <sundials/sundials_types.h>     // definition of realtype
#include <sundials/sundials_dense.h>     // generic DENSE solver
#include <sundials/sundials_types.h>     // definition of realtype
#include <sunlinsol/sunlinsol_dense.h>   // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_dense.h>   // access to dense SUNMatrix


// Declare and define file-level constants
const double T = 625 + 273;
const double G = 3e-3;
const double he_2_dpa = 5e-6;
// double k_B;
// double gamma_b;
// double Omega;
// double Ef_v;
// double B;

std::map<std::string, double> props;

std::map<std::string, double> get_properties(char *argv[]) {
// std::map<std::string, double> get_properties() {
    // std::map<std::string, double> props;
    // props["nu_v"] = std::stod(argv[1]);
    // props["nu_i"] = std::stod(argv[2]);
    // props["nu_g"] = std::stod(argv[3]);
    // props["Em_v"] = std::stod(argv[4]);
    // props["Em_i"] = std::stod(argv[5]);
    // props["Em_g"] = std::stod(argv[6]);
    // props["Eb_v_g"] = std::stod(argv[7]);
    // props["Eb_v_2g"] = std::stod(argv[8]);
    // props["Eb_2g"] = std::stod(argv[9]);
    // props["Ef_v"] = std::stod(argv[10]);
    // props["a0"] = std::stod(argv[11]);
    // props["Omega"] = std::stod(argv[12]);
    // props["f"] = std::stod(argv[13]);
    // props["b"] = std::stod(argv[14]);
    // props["k_B"] = std::stod(argv[15]);
    // props["gamma_b"] = std::stod(argv[16]);
    // props["B"] = std::stod(argv[17]);
    // props["r_ppt"] = std::stod(argv[18]);
    // props["d"] = std::stod(argv[19]);
    // props["N_ppt"] = std::stod(argv[20]);
    // props["rho"] = std::stod(argv[21]);
    // props["Z_i"] = std::stod(argv[22]);


    std::string nu_v_str = argv[1];
    std::string nu_i_str = argv[2];
    std::string nu_g_str = argv[3];

    std::string Em_v_str = argv[4];
    std::string Em_i_str = argv[5];
    std::string Em_g_str = argv[6];
    std::string Eb_v_g_str = argv[7];
    std::string Eb_v_2g_str = argv[8];
    std::string Eb_2g_str = argv[9];
    std::string Ef_v_str = argv[10];
    std::string a0_str = argv[11];
    std::string Omega_str = argv[12];
    std::string f_str = argv[13];
    std::string b_str = argv[14];
    std::string k_B_str = argv[15];
    std::string gamma_b_str = argv[16];
    std::string B_str = argv[17];
    std::string r_ppt_str = argv[18];
    std::string d_str = argv[19];
    std::string N_ppt_str = argv[20];
    std::string rho_str = argv[21];
    std::string Z_i_str = argv[22];

    double nu_v = std::stod(nu_v_str);
    double nu_i = std::stod(nu_i_str);
    double nu_g = std::stod(nu_g_str);
    double Em_v = std::stod(Em_v_str);
    double Em_i = std::stod(Em_i_str);
    double Em_g = std::stod(Em_g_str);
    double Eb_v_g = std::stod(Eb_v_g_str);
    double Eb_v_2g = std::stod(Eb_v_2g_str);
    double Eb_2g = std::stod(Eb_2g_str);
    double Ef_v = std::stod(Ef_v_str);
    double a0 = std::stod(a0_str);
    double Omega = std::stod(Omega_str);
    double f = std::stod(f_str);
    double b = std::stod(b_str);
    double k_B = std::stod(k_B_str);
    double gamma_b = std::stod(gamma_b_str);
    double B = std::stod(B_str);
    double r_ppt = std::stod(r_ppt_str);
    double d = std::stod(d_str);
    double N_ppt = std::stod(N_ppt_str);
    double rho = std::stod(rho_str);
    double Z_i = std::stod(Z_i_str);

    // double nu_v = 5000000000000.0;
    // double nu_i = 50000000000000.0;
    // double nu_g = 50000000000000.0;
    // double Em_v = 1.4;
    // double Em_i = 0.2;
    // double Em_g = 0.2;
    // double Eb_v_g = 2.4;
    // double Eb_v_2g = 3.5;
    // double Eb_2g = 0.79;
    // double Ef_v = 1.6;
    // double a0 = 3.63e-10;
    // double Omega = 1.1958036749999998e-29;
    // double f = 0.1;
    // double b = 10.0;
    // double k_B = 8.617e-05;
    // double gamma_b = 6.24e+18;
    // double B = 1.75e-29;
    // double r_ppt = 1e-08;
    // double d = 3e-05;
    // double N_ppt = 1e+16;
    // double rho = 300000000000000.0;
    // double Z_i = 1.2;

    // Reaction frequencies
    double alpha = 48 * nu_i * exp(-Em_i / (k_B * T));
    double beta = 48 * nu_g * exp(-Em_g / (k_B * T));
    double gamma = 48 * nu_v * exp(-Em_v / (k_B * T));

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
    double P = f * G;
    double delta = b * G;
    double G_He = he_2_dpa * G;

    return {
        {"alpha", alpha}, {"beta", beta}, {"gamma", gamma}, {"e1", e1},
        {"e2", e2}, {"e4", e4}, {"e5", e5}, {"delta", delta}, {"D_i", D_i},
        {"D_v", D_v}, {"D_g", D_g}, {"P", P}, {"G_He", G_He}, {"C_v_e", C_v_e},
        {"Omega", Omega}, {"C_ppt", C_ppt}, {"rho", rho}, {"a0", a0}, {"d", d},
        {"Z_i", Z_i}, {"r_ppt", r_ppt}, {"k_B", k_B}, {"gamma_b", gamma_b}, {"B", B},
        {"T", T}, {"Ef_v", Ef_v} 
    };
}

std::tuple<double, double, double> get_bubble_props(double R, double m, std::map<std::string, double>& props) {
    // double k_B = p_value["k_B"];
    // double T = ::T;
    // double B = p_value["B"];
    // double Omega = p_value["Omega"];
    // double Ef_v = p_value["Ef_v"];
    // double gamma_b = p_value["gamma_b"];
    // double a0 = p_value["a0"];

    double pressure = m * props["k_B"] * props["T"] / (4.0 * M_PI * R * R * R / 3.0 - m * props["B"]);
    double work = (2 * props["gamma_b"] / R - pressure) * props["Omega"];
    double Eb_v_B = props["Ef_v"] + work;
    double e3 = exp(-Eb_v_B / (props["k_B"] * props["T"]));
    double epsilon = (4 * M_PI / 48) * (R / props["a0"]);

    return std::make_tuple(e3, epsilon, work);
}

// Right-hand side function
int rhs(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data) {

    // Debug print to check if user_data is null
    if (user_data == nullptr) {
        std::cerr << "Error: user_data is null." << std::endl;
        return -1;
    }

    // Unpack the state vector
    // std::vector<double> y_vec(NV_LENGTH_S(y));
    // for (int i = 0; i < NV_LENGTH_S(y); ++i) {
    //     y_vec[i] = NV_Ith_S(y, i);
    // }

    Parameters* params = static_cast<Parameters*>(user_data);
    std::map<std::string, double>& props = params->props;
    // double T = params-> T;
    // // 625 + 273;
    // double G = 3e-3;
    // double he_2_dpa = 5e-6;
    // double k_B;
    // double gamma_b;
    // double Omega;
    // double Ef_v;
    // double B;
    // Unpack parameters (props)
    double alpha = props["alpha"];
    double beta = props["beta"];
    double gamma_val = props["gamma"];
    double e1 = props["e1"];
    double e2 = props["e2"];
    double e4 = props["e4"];
    double e5 = props["e5"];
    double delta = props["delta"];
    double D_i = props["D_i"];
    double D_v = props["D_v"];
    double D_g = props["D_g"];
    double P = props["P"];
    double G_He = props["G_He"];
    double C_v_e = props["C_v_e"];
    double Omega = props["Omega"];
    double C_ppt = props["C_ppt"];
    double rho = props["rho"];
    double a0 = props["a0"];
    double d = props["d"];
    double Z_i = props["Z_i"];
    double r_ppt = props["r_ppt"];

    // Prepare time-dependent bubble-related parameters
    auto [e3, epsilon, work] = get_bubble_props(NV_Ith_S(y, 9), NV_Ith_S(y, 8), props);
    auto [e3_prime, epsilon_ppt, work_ppt] = get_bubble_props(NV_Ith_S(y, 11), NV_Ith_S(y, 10), props);

    double k_square = 4 * M_PI * NV_Ith_S(y, 9) * NV_Ith_S(y, 7) / Omega + rho;
    double C_gb = a0 * a0 * sqrt(k_square) / (8 * d);

    // Sink concentrations
    double bubble_sink = 4 * M_PI * NV_Ith_S(y, 9) * NV_Ith_S(y, 7) / Omega;
    double precipitate_sink = 4 * M_PI * NV_Ith_S(y, 11) * C_ppt / Omega;

    double C_s_v = (a0 * a0 / 48) * (rho + bubble_sink);
    double C_s_i = (a0 * a0 / 48) * (Z_i * rho + bubble_sink);

    // Compute derivatives directly using NV_Ith_S(y, i)
    NV_Ith_S(ydot, 0) = (P + (beta * e1 + delta) * NV_Ith_S(y, 3) - (alpha * NV_Ith_S(y, 1) + beta * NV_Ith_S(y, 2) + gamma_val * (C_s_v + NV_Ith_S(y, 3) + 2 * (NV_Ith_S(y, 5) + NV_Ith_S(y, 4)) + 3 * NV_Ith_S(y, 6))) * NV_Ith_S(y, 0));
    NV_Ith_S(ydot, 1) = P - alpha * (NV_Ith_S(y, 0) + NV_Ith_S(y, 3) + 2 * NV_Ith_S(y, 4) + 3 * NV_Ith_S(y, 6) + C_s_i) * NV_Ith_S(y, 1);
    NV_Ith_S(ydot, 2) = (G_He - beta * NV_Ith_S(y, 2) * (NV_Ith_S(y, 0) + 2 * NV_Ith_S(y, 2) + NV_Ith_S(y, 3) + 2 * NV_Ith_S(y, 5) + 2 * NV_Ith_S(y, 4) + C_gb + epsilon * NV_Ith_S(y, 7)) + delta * (NV_Ith_S(y, 3) + 2 * NV_Ith_S(y, 4) + 4 * NV_Ith_S(y, 5) + 3 * NV_Ith_S(y, 6) + NV_Ith_S(y, 8) * NV_Ith_S(y, 7) + NV_Ith_S(y, 12) + NV_Ith_S(y, 10) * C_ppt) + alpha * NV_Ith_S(y, 1) * NV_Ith_S(y, 3) + beta * (e1 * NV_Ith_S(y, 3) + e2 * NV_Ith_S(y, 4)));
    NV_Ith_S(ydot, 3) = (beta * NV_Ith_S(y, 2) * NV_Ith_S(y, 0) + (beta * e2 + 2 * delta) * NV_Ith_S(y, 4) - (beta * e1 + delta + alpha * NV_Ith_S(y, 1) + beta * NV_Ith_S(y, 2)) * NV_Ith_S(y, 3));
    NV_Ith_S(ydot, 4) = (beta * NV_Ith_S(y, 2) * NV_Ith_S(y, 3) + 3 * delta * NV_Ith_S(y, 6) + 2 * gamma_val * NV_Ith_S(y, 0) * NV_Ith_S(y, 5) - (2 * beta * NV_Ith_S(y, 2) + 2 * delta + beta * e2 + 2 * alpha * NV_Ith_S(y, 1)) * NV_Ith_S(y, 4));
    NV_Ith_S(ydot, 5) = (alpha * NV_Ith_S(y, 1) * NV_Ith_S(y, 4) + 2 * beta * NV_Ith_S(y, 2) * NV_Ith_S(y, 2) - (2 * delta + 2 * gamma_val * NV_Ith_S(y, 0) + 2 * beta * NV_Ith_S(y, 2)) * NV_Ith_S(y, 5));
    NV_Ith_S(ydot, 6) = (2 * beta * NV_Ith_S(y, 2) * (NV_Ith_S(y, 4) + NV_Ith_S(y, 5)) - 3 * NV_Ith_S(y, 6) * (delta + alpha * NV_Ith_S(y, 1) + beta * NV_Ith_S(y, 2) + gamma_val * NV_Ith_S(y, 0)));
    NV_Ith_S(ydot, 7) = ((12 * beta * NV_Ith_S(y, 2) + 9 * gamma_val * NV_Ith_S(y, 0)) * NV_Ith_S(y, 6)) / NV_Ith_S(y, 8);
    NV_Ith_S(ydot, 8) = epsilon * beta * NV_Ith_S(y, 2) - delta * NV_Ith_S(y, 8);
    NV_Ith_S(ydot, 9) = (a0 * a0 / NV_Ith_S(y, 9)) * (gamma_val * NV_Ith_S(y, 0) - alpha * NV_Ith_S(y, 1) - gamma_val * (e3 - e4));
    if (NV_Ith_S(ydot, 9) < 0) NV_Ith_S(ydot, 9) = 0;
    NV_Ith_S(ydot, 10) = epsilon_ppt * beta * NV_Ith_S(y, 2) - delta * NV_Ith_S(y, 10);
    double r_p_equiv = sqrt(NV_Ith_S(y, 11) * NV_Ith_S(y, 11) + r_ppt * r_ppt);
    NV_Ith_S(ydot, 11) = (a0 * a0 / r_p_equiv) * (gamma_val * NV_Ith_S(y, 0) - alpha * NV_Ith_S(y, 1) - gamma_val * (e3_prime - e4));
    if (NV_Ith_S(ydot, 11) < 0) NV_Ith_S(ydot, 11) = 0;
    NV_Ith_S(ydot, 12) = beta * C_gb * NV_Ith_S(y, 2) - delta * NV_Ith_S(y, 12);

    // Print values to console
    // std::cout << "t: " << t << std::endl;
    // std::cout << "ydot: ";
    // for (int i = 0; i < NV_LENGTH_S(ydot); ++i) {
    //     std::cout << NV_Ith_S(ydot, i) << " ";
    // }
    // std::cout << std::endl;
    return 0;
}

int main(int argc, char *argv[]) {

    auto props = get_properties(argv);
    // auto props = get_properties();

    Parameters params;
    params.props = props;

    // for (auto const& [key, val] : props) {
    //     std::cout << key << ": " << val << std::endl;
    // }

    // Problem setup
    sunrealtype t0 = 1e-6;
    sunrealtype tf = 1e6;
    sunrealtype reltol = 1e-8;
    sunrealtype abstol = 1e-20;
    sunindextype neq = 13; // number of equations
    // Cast sunindextype to int for loop usage
    int neq_int = static_cast<int>(neq);

    

    // Create SUNContext
    SUNContext sunctx;
    SUNContext_Create(SUN_COMM_NULL, &sunctx);

    // Define the initial state vector 'y'
    double floor_val = 1e-20;
    N_Vector y = N_VNew_Serial(neq, sunctx);
    for (int i = 0; i < neq; ++i) NV_Ith_S(y, i) = floor_val;
    NV_Ith_S(y, 8) = 2.0;
    NV_Ith_S(y, 9) = 5e-10;
    NV_Ith_S(y, 10) = 2.0;
    NV_Ith_S(y, 11) = 5e-10;

    // Create CVODE solver using BDF method for stiff problems
    void *cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (cvode_mem == NULL) {
        std::cerr << "Error in CVodeCreate" << std::endl;
        return 1;
    }
    // std::cout << "CVODE solver created successfully" << std::endl;

    // Set user data
    int flag = CVodeSetUserData(cvode_mem, &params);
    if (flag != CV_SUCCESS) {
        std::cerr << "Error in CVodeSetUserData" << std::endl;
        return 1;
    }
    // std::cout << "User data set successfully" << std::endl;

    // Initialize CVODE solver
    flag = CVodeInit(cvode_mem, rhs, t0, y);
    if (flag != CV_SUCCESS) {
        std::cerr << "Error in CVodeInit" << std::endl;
        return 1;
    }
    // std::cout << "CVODE solver initialized successfully" << std::endl;

    // Set max steps
    flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
    if (flag != CV_SUCCESS) {
        std::cerr << "Error in CVodeSetMaxNumSteps" << std::endl;
        return 1;
    }

    // Set Stability Limit Detection
    flag = CVodeSetStabLimDet(cvode_mem, SUNTRUE);
    if (flag != CV_SUCCESS) {
        std::cerr << "Error in CVodeSetStabLimDet" << std::endl;
        return 1;
    }

    // Set minimum step size
    flag = CVodeSetMinStep(cvode_mem, 1e-20);
    if (flag != CV_SUCCESS) {
        std::cerr << "Error in CVodeSetMinStep" << std::endl;
        return 1;
    }
    // Set Initial Step Size
    // flag = CVodeSetInitStep(cvode_mem, 1e-12);
    // if (flag != CV_SUCCESS) {
    //     std::cerr << "Error in CVodeSetInitStep" << std::endl;
    //     return 1;
    // }
   // Set tolerances
    flag = CVodeSStolerances(cvode_mem, reltol, abstol);
    if (flag != CV_SUCCESS) {
        std::cerr << "Error in CVodeSStolerances" << std::endl;
        return 1;
    }
    // std::cout << "Tolerances set successfully" << std::endl;

    // Create dense matrix for the Jacobian
    SUNMatrix A = SUNDenseMatrix(neq, neq, sunctx);
    if (A == NULL) {
        std::cerr << "Error creating SUNDenseMatrix" << std::endl;
        return 1;
    }

    // // Create dense linear solver object
    SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);
    if (LS == NULL) {
        std::cerr << "Error creating SUNLinSol_Dense" << std::endl;
        return 1;
    }

    // // Attach the matrix and linear solver to CVODE
    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (flag != CV_SUCCESS) {
        std::cerr << "Error in CVodeSetLinearSolver" << std::endl;
        return 1;
    }
     

    // Define time points
    double t_span_start = 1e-6;
    double t_span_end = 1e6;
    int time_points = 200;

    std::vector<double> t_eval(time_points);
    double log_t_span_start = std::log10(t_span_start);
    double log_t_span_end = std::log10(t_span_end);
    double step = (log_t_span_end - log_t_span_start) / (t_eval.size() - 1);
    for (size_t i = 0; i < t_eval.size(); ++i) {
        t_eval[i] = std::pow(10, log_t_span_start + i * step);
    }

    sunrealtype t = t0;

    // Main time-stepping loop
    std::vector<std::vector<double>> results;
    // Push initial state
    results.push_back({NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2), NV_Ith_S(y, 3), NV_Ith_S(y, 4), NV_Ith_S(y, 5), NV_Ith_S(y, 6), NV_Ith_S(y, 7), NV_Ith_S(y, 8), NV_Ith_S(y, 9), NV_Ith_S(y, 10), NV_Ith_S(y, 11), NV_Ith_S(y, 12)});
    
    for (size_t i = 1; i < t_eval.size(); ++i) {
        sunrealtype tout = t_eval[i];
        flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        if (flag < 0) {
            std::cerr << "Error in CVode at time " << t << " with flag " << flag << std::endl;
            break;
        }
        std::vector<double> state(neq_int);
        for (int j = 0; j < neq_int; ++j) {
            state[j] = NV_Ith_S(y, j);
        }
        results.push_back(state);
    }

    // std::cout << "Integration complete" << std::endl;

    // Set output precision to scientific notation
    std::cout << std::scientific << std::setprecision(10);
    // Output results as a 200 x 13 matrix
    for (const auto& row : results) {
        for (const auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }

    // Free memory
    N_VDestroy(y);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    CVodeFree(&cvode_mem);
    SUNContext_Free(&sunctx);

    return 0;
}