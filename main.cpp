#include "parameters.h"
#include <iomanip>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES
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

#include <chrono>

// Declare and define file-level constants
const double T = 550 + 273;
// const double G = 3e-3;
// const double he_2_dpa = 5e-6;
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

    // double tau_oct = 150000000.0;
    // double b = 2e-10;
    // double D_o = 2e-05;
    // double SFE = 0.2;
    // double delta = 4e-09;
    // double sigma_o = 20000000.0;
    // double E = 177000000000.0;
    // double nu = 0.3;
    // double N_p = 1.63e+16;
    // double r_p = 4.05e-08;
    // double k = 1.38065e-23;
    // double Omega = 1.19e-29;
    // double W_g = 6e-19;
    // double a1 = 100000000000.0;
    // double c_jog = 0.0389;
    // double E_core = 2.0826e-19;
    // double E_s = 4.4856e-19;
    // double E_m = 2.2428e-19;
    // double K_c = 10.0;
    // double Beta = 115000.0;
    // double xi = 1.0;
    // double zeta = 0.425;

    std::string tau_oct_str = argv[1];
    std::string b_str = argv[2];
    std::string D_o_str = argv[3];
    std::string SFE_str = argv[4];
    std::string delta_str = argv[5];
    std::string sigma_o_str = argv[6];
    std::string E_str = argv[7];
    std::string nu_str = argv[8];
    std::string N_p_str = argv[9];
    std::string r_p_str = argv[10];
    std::string k_str = argv[11];
    std::string Omega_str = argv[12];
    std::string W_g_str = argv[13];
    std::string a1_str = argv[14];
    std::string c_jog_str = argv[15];
    std::string E_core_str = argv[16];
    std::string E_s_str = argv[17];
    std::string E_m_str = argv[18];
    std::string K_c_str = argv[19];
    std::string Beta_str = argv[20];
    std::string xi_str = argv[21];
    std::string zeta_str = argv[22];

    double tau_oct = std::stod(tau_oct_str);
    double b = std::stod(b_str);
    double D_o = std::stod(D_o_str);
    double SFE = std::stod(SFE_str);
    double delta = std::stod(delta_str);
    double sigma_o = std::stod(sigma_o_str);
    double E = std::stod(E_str);
    double nu = std::stod(nu_str);
    double N_p = std::stod(N_p_str);
    double r_p = std::stod(r_p_str);
    double k = std::stod(k_str);
    double Omega = std::stod(Omega_str);
    double W_g = std::stod(W_g_str);
    double a1 = std::stod(a1_str);
    double c_jog = std::stod(c_jog_str);
    double E_core = std::stod(E_core_str);
    double E_s = std::stod(E_s_str);
    double E_m = std::stod(E_m_str);
    double K_c = std::stod(K_c_str);
    double Beta = std::stod(Beta_str);
    double xi = std::stod(xi_str);
    double zeta = std::stod(zeta_str);

    double mu = E / (2 * (1 + nu));
    double eta_v = 1e3 * 1 * c_jog * b * (SFE / (mu * b)) * (SFE / (mu * b));
    double D_p = D_o * exp(-E_core / (k * T));
    double D_v = D_o * exp(-E_m / (k * T));
    double D_s = D_o * exp(-E_s / (k * T));

    return {
        {"tau_oct", tau_oct}, {"b", b}, {"D_o", D_o}, {"SFE", SFE}, {"delta", delta},
        {"sigma_o", sigma_o}, {"E", E}, {"nu", nu}, {"N_p", N_p}, {"r_p", r_p},
        {"k", k}, {"Omega", Omega}, {"W_g", W_g}, {"a1", a1}, {"c_jog", c_jog},
        {"E_core", E_core}, {"E_s", E_s}, {"E_m", E_m}, {"K_c", K_c}, {"Beta", Beta},
        {"xi", xi}, {"zeta", zeta}, {"mu", mu}, {"eta_v", eta_v}, {"D_p", D_p},
        {"D_v", D_v}, {"D_s", D_s}
    };
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

    double tau_oct = props["tau_oct"];
    double b = props["b"];
    double D_o = props["D_o"];
    double SFE = props["SFE"];
    double delta = props["delta"];
    double sigma_o = props["sigma_o"];
    double E = props["E"];
    double nu = props["nu"];
    double N_p = props["N_p"];
    double r_p = props["r_p"];
    double k = props["k"];
    double Omega = props["Omega"];
    double W_g = props["W_g"];
    double a1 = props["a1"];
    double c_jog = props["c_jog"];
    double E_core = props["E_core"];
    double E_s = props["E_s"];
    double E_m = props["E_m"];
    double K_c = props["K_c"];
    double Beta = props["Beta"];
    double xi = props["xi"];
    double zeta = props["zeta"];
    double mu = props["mu"];
    double eta_v = props["eta_v"];
    double D_p = props["D_p"];
    double D_v = props["D_v"];
    double D_s = props["D_s"];

    double rho_m = NV_Ith_S(y, 0);
    double rho_s = NV_Ith_S(y, 1);
    double rho_b = NV_Ith_S(y, 2);
    double R_sb = NV_Ith_S(y, 3);

    double h = 1/(R_sb*(rho_s + rho_b)); // [m] dislocations spacing within the subgrain walls(eq3 on 1990-ghon-matt paper)
    double gamma_sb = mu*b*b*rho_b*R_sb/3; // [J/m^2] Low-angle subgrain boundry energy per unit area(eq36 on 1990-ghon-matt paper)
    double p_s = (4.0/3.0)*mu*pow(b, 2)*rho_b; // [Pa] Pressure for subgrain growth(Gibbs-Thompson effect)(eq38 on 1990-ghon-matt paper)
    double M_sb_c = 2*M_PI*b*D_p*Omega/(h*h*k*T); // unit?[m^3/(N.s)] Core mobility(eq39 on 1990-ghon-matt paper)
    double M_sb_L = 2*M_PI*eta_v*D_v*Omega/(b*k*T); // unit?[m^3/(N.s)] Lattice mobility(eq40 on 1990-ghon-matt paper)
    double M_sb = M_sb_c + M_sb_L; // unir?[m^3/(N.s)]
    double lambda_d = 1/sqrt(rho_m); // [m] inter-dislocation spacing
    double lambda_p = 1/sqrt(N_p*r_p); // [m] inter-precipitates spacing
    double lambda = 1/(1/lambda_d + 1/lambda_p); // [m] inter-obstacles spacing considering both dislocations and precipitates
    double alpha = 1/(M_PI*(1-nu));
    double sigma_i = mu*b/(2*M_PI*lambda) + xi*mu*b*sqrt(rho_s); // [Pa] long range internal stress(eq7 on 1990-ghon-matt paper)
    double sigma_e = tau_oct - sigma_i - sigma_o; // [Pa] Effective stress on dislocations
    double sigma_s = mu*b*sqrt(rho_s)*alpha; // [Pa] Static Dislocations stress???
    double sigma_m = mu*b*sqrt(rho_m)*alpha; // [Pa] Mobile Dislocations stress???
    double v_g = (sigma_e*Omega*a1/(k*T))*exp(-W_g/(k*T)); // [?] Empirical glide velocity (eq12 on 1990-ghon-matt paper)
    double v_cs = c_jog*2*M_PI*(D_s/b)*(Omega*sigma_s/(k*T))/log(1/sqrt(b/R_sb));
    double v_cm = c_jog*2*M_PI*(D_s/b)*(Omega*sigma_m/(k*T))/log(1/sqrt(b*b*rho_m));

    double R1 = pow(sqrt(rho_m), 3);
    double R2 = Beta/(h*h*R_sb);
    double R3 = rho_m/(2*R_sb);
    double R4 = 8*pow(sqrt(rho_m), 3)*v_cm/v_g;
    double R5 = delta*rho_m*(rho_m + rho_s);
    double R6 = 8.0*(rho_s/h)*(v_cs/v_g);
    double R7 = delta*rho_m*rho_s;
    double R8 = 8.0*(1.0 - 2*zeta)*rho_s*(v_cs/h);
    double R9 = (rho_b/R_sb)*M_sb*(p_s - 2*M_PI*pow(r_p, 2)*N_p*gamma_sb);
    double R10 = M_sb*(p_s - 2*M_PI*r_p*r_p*N_p*gamma_sb);
    double R11 = mu*eta_v*K_c*R_sb*(sqrt(rho_m + rho_s) - K_c/(2*R_sb))*(Omega*D_s/(k*T));

    // Compute derivaties

    NV_Ith_S(ydot, 0) = v_g*(R1 + R2 - R3 - R4 - R5); // rho_m_dot
    NV_Ith_S(ydot, 1) = v_g*(R3 - R6 - R7); // rho_s_dot
    NV_Ith_S(ydot, 2) = R8 - R9; // rho_b_dot
    NV_Ith_S(ydot, 3) = R10 - R11; // R_sb_dot
    NV_Ith_S(ydot, 4) = rho_m*b*v_g;

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
    sunrealtype t0 = 1e-10;
    sunrealtype tf = 4e7;
    sunrealtype reltol = 1e-3;
    sunrealtype abstol = 1e-6;
    sunindextype neq = 5; // number of equations
    // Cast sunindextype to int for loop usage
    int neq_int = static_cast<int>(neq);

    

    // Create SUNContext
    SUNContext sunctx;
    SUNContext_Create(SUN_COMM_NULL, &sunctx);

    // Define the initial state vector 'y'
    N_Vector y = N_VNew_Serial(neq, sunctx);
    
    NV_Ith_S(y, 0) = 1e14;
    NV_Ith_S(y, 1) = 1e11;
    NV_Ith_S(y, 2) = 1e11;
    NV_Ith_S(y, 3) = 10e-6;
    NV_Ith_S(y, 4) = 0.0;

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
    // flag = CVodeSetStabLimDet(cvode_mem, SUNTRUE);
    // if (flag != CV_SUCCESS) {
    //     std::cerr << "Error in CVodeSetStabLimDet" << std::endl;
    //     return 1;
    // }

    // Set minimum step size
    // flag = CVodeSetMinStep(cvode_mem, 1e-20);
    // if (flag != CV_SUCCESS) {
    //     std::cerr << "Error in CVodeSetMinStep" << std::endl;
    //     return 1;
    // }
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
    double t_span_start = 1e-10;
    double t_span_end = 4e7;
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
    results.push_back({NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2), NV_Ith_S(y, 3), NV_Ith_S(y, 4)});
    
    // using namespace std::chrono;
    // auto start = std::chrono::high_resolution_clock::now();

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

    // auto end = std::chrono::high_resolution_clock::now();
    // // Calculate the duration in milliseconds
    // std::chrono::duration<double, std::milli> duration = end - start;
    // // Output the duration
    // std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;

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