// parameters.h
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <map>
#include <string>

// Define a structure to hold all necessary parameters
struct Parameters {
    // double T = 625 + 273;
    // double G = 3e-3;
    // double he_2_dpa = 5e-6;
    // double k_B;
    // double gamma_b;
    // double Omega;
    // double Ef_v;
    // double B;
    std::map<std::string, double> props;
};

#endif // PARAMETERS_H