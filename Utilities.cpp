#include "Utilities.h"
#include <OpenXLSX/OpenXLSX.hpp>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <tuple>

// Define global variables
// state_type y;
// double floor_val = 1e-20;
// int array_length = 13;
// vector_type y(array_length, floor_val);


// Define file-scope variables
static int font1 = 16;
static int font2 = 18;
static double T = 625 + 273;
static double G = 3e-3;
static double he_2_dpa = 5e-6;

static std::map<std::string, double> p_value;

// // Define state type
// typedef std::vector<double> state_type;

static double nu_v;
static double nu_i;
static double nu_g;
static double Em_v;
static double Em_i;
static double Em_g;
static double Eb_v_g;
static double Eb_v_2g;
static double Eb_2g;
static double Ef_v;
static double a0;
static double Omega;
static double f;
static double b;
static double k_B;
static double gamma_b;
static double B;
static double r_ppt;
static double d;
static double N_ppt;
static double rho;
static double Z_i;

void read_excel(const std::string& file_path) {
    // Open the Excel file
    OpenXLSX::XLDocument doc;
    doc.open(file_path);
    auto wks = doc.workbook().worksheet("steel");

    // Read the values from the Excel file
    for (auto& row : wks.rows(22)) {
        std::vector<OpenXLSX::XLCellValue> row_values = row.values();
        std::string variable_name = row_values[0].get<std::string>();
        int index_val = 0;

        for (auto& value : row_values) {
            
            if (index_val == 1 ) {
                if (value.type() == OpenXLSX::XLValueType::Integer) {
                    // Cast the value as double and save to p_value
                    int64_t int_value = value.get<int64_t>();
                    p_value[variable_name] = static_cast<double>(int_value);
                }
                // else if (value.type() == OpenXLSX::XLValueType::Float) {
                //     // Cast the value as double and save to p_value
                //     p_value[variable_name] = static_cast<double>(value.get<float>());
                // }
                else {
                    p_value[variable_name] = value.get<double>();
                }
                break;
            }
            index_val++;
        }
        // p_value[variable_name] = variable_value;
    }
    doc.close();

    nu_v = p_value["nu_v"];
    nu_i = p_value["nu_i"];
    nu_g = p_value["nu_g"];
    Em_v = p_value["Em_v"];
    Em_i = p_value["Em_i"];
    Em_g = p_value["Em_g"];
    Eb_v_g = p_value["Eb_v_g"];
    Eb_v_2g = p_value["Eb_v_2g"];
    Eb_2g = p_value["Eb_2g"];
    Ef_v = p_value["Ef_v"];
    a0 = p_value["a0"];
    Omega = p_value["Omega"];
    f = p_value["f"];
    b = p_value["b"];
    k_B = p_value["k_B"];
    gamma_b = p_value["gamma_b"];
    B = p_value["B"];
    r_ppt = p_value["r_ppt"];
    d = p_value["d"];
    N_ppt = p_value["N_ppt"];
    rho = p_value["rho"];
    Z_i = p_value["Z_i"];

}

std::map<std::string, double> get_props() {
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
        {"Z_i", Z_i}, {"r_ppt", r_ppt}
    };
}

std::tuple<double, double, double> get_bubble_props(double R, double m) {
    // double k_B = p_value["k_B"];
    // double T = ::T;
    // double B = p_value["B"];
    // double Omega = p_value["Omega"];
    // double Ef_v = p_value["Ef_v"];
    // double gamma_b = p_value["gamma_b"];
    // double a0 = p_value["a0"];

    double pressure = m * k_B * T / (4.0 * M_PI * R * R * R / 3.0 - m * B);
    double work = (2 * gamma_b / R - pressure) * Omega;
    double Eb_v_B = Ef_v + work;
    double e3 = exp(-Eb_v_B / (k_B * T));
    double epsilon = (4 * M_PI / 48) * (R / a0);

    return std::make_tuple(e3, epsilon, work);
}

void initialize() {
    // floor_val = 1e-20;
    // array_length = 13;
    
    // y.resize(array_length, floor_val);
    y[8] = 2;
    y[9] = 5e-10;
    y[10] = 2;
    y[11] = 5e-10;
    // return std::make_tuple(y0_vec, floor_val, array_length);
}