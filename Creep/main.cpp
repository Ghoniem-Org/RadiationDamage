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

std::map<std::string, double> get_properties(int argc, char *argv[]) {
    std::map<std::string, double> props;

    // Loop through all arguments
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
    double t0 = props["t0"];
    double tf = props["tf"];
    double time_points = props["time_points"];

    double mu = E / (2 * (1 + nu));
    double eta_v = 1e3 * 1 * c_jog * b * (SFE / (mu * b)) * (SFE / (mu * b));
    double D_p = D_o * exp(-E_core / (k * T));
    double D_v = D_o * exp(-E_m / (k * T));
    double D_s = D_o * exp(-E_s / (k * T));
    
    props["mu"] = mu;
    props["eta_v"] = eta_v;
    props["D_p"] = D_p;
    props["D_v"] = D_v;
    props["D_s"] = D_s;

    return props;
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

    auto props = get_properties(argc, argv);
    // auto props = get_properties();

    Parameters params;
    params.props = props;

    // for (auto const& [key, val] : props) {
    //     std::cout << key << ": " << val << std::endl;
    // }

    // Problem setup
    sunrealtype t0 = props["t0"];
    sunrealtype tf = props["tf"];
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
    double t_span_start = props["t0"];
    double t_span_end = props["tf"];
    int time_points = int(props["time_points"]);

    // Define time vector in linear scale
    std::vector<double> t_eval(time_points);
    double step = (t_span_end - t_span_start) / (time_points - 1);
    for (size_t i = 0; i < t_eval.size(); ++i) {
        t_eval[i] = t_span_start + i * step;
    }
    // std::vector<double> t_eval(time_points);
    // double log_t_span_start = std::log10(t_span_start);
    // double log_t_span_end = std::log10(t_span_end);
    // double step = (log_t_span_end - log_t_span_start) / (t_eval.size() - 1);
    // for (size_t i = 0; i < t_eval.size(); ++i) {
    //     t_eval[i] = std::pow(10, log_t_span_start + i * step);
    // }

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