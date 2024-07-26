#include "Derivatives.h"
#include "Utilities.h"
#include <vector>
#include <cmath>
#include <map>
#include <iostream>

// // Define global variables
std::map<std::string, double> props;

// // Define state type
// typedef std::vector<double> state_type;

// File-scope variables
static double P;
static double alpha;
static double beta;
static double gamma_val;
static double delta;
static double e1;
static double e2;
static double e4;
static double e5;
static double Omega;
static double G_He;
static double C_ppt;
static double rho;
static double a0;
static double d;
static double Z_i;
static double r_ppt;

static double e3, epsilon, work;
static double e3_prime, epsilon_ppt, work_ppt;
static double k_square;
static double C_gb;
static double bubble_sink;
static double precipitate_sink;
static double C_s_v;
static double C_s_i;

template< class State >
void load_variables(const State &y) {
    props = get_props();

    P = props["P"];
    alpha = props["alpha"];
    beta = props["beta"];
    gamma_val = props["gamma"];
    delta = props["delta"];
    e1 = props["e1"];
    e2 = props["e2"];
    e4 = props["e4"];
    e5 = props["e5"];
    Omega = props["Omega"];
    G_He = props["G_He"];
    C_ppt = props["C_ppt"];
    rho = props["rho"];
    a0 = props["a0"];
    d = props["d"];
    Z_i = props["Z_i"];
    r_ppt = props["r_ppt"];

    // int array_length = y.size();
    // dydt.resize(array_length);
    // const double &C_v = y[0];
    // const double &C_i = y[1];
    // const double &C_g = y[2];
    // const double &C_gv = y[3];
    // const double &C_2gv = y[4];
    // const double &C_2g = y[5];
    // const double &C_star = y[6];
    // const double &C_b = y[7];
    // const double &m = y[8];
    // const double &R = y[9];
    // const double &m_ppt = y[10];
    // const double &R_pt = y[11];
    // const double &M_gb = y[12];

    // Prepare time-dependent bubble-related parameters:
    std::tie(e3, epsilon, work) = get_bubble_props(y[9], y[8]);

    std::tie(e3_prime, epsilon_ppt, work_ppt) = get_bubble_props(y[11], y[10]);

    k_square = 4 * M_PI * y[9] * y[7] / Omega + rho;
    C_gb = a0 * a0 * sqrt(k_square) / (8 * d);

    // Sink concentrations
    bubble_sink = 4 * M_PI * y[9] * y[7] / Omega;
    precipitate_sink = 4 * M_PI * y[11] * C_ppt / Omega;

    C_s_v = (a0 * a0 / 48) * (rho + bubble_sink);
    C_s_i = (a0 * a0 / 48) * (Z_i * rho + bubble_sink);

}

// typedef boost::numeric::ublas::vector< double > vector_type;
// typedef boost::numeric::ublas::matrix< double > matrix_type;
template<class State>
void derivatives::operator()( const State &y , State &dydt , double t )
    {   
        props = get_props();

        P = props["P"];
        alpha = props["alpha"];
        beta = props["beta"];
        gamma_val = props["gamma"];
        delta = props["delta"];
        e1 = props["e1"];
        e2 = props["e2"];
        e4 = props["e4"];
        e5 = props["e5"];
        Omega = props["Omega"];
        G_He = props["G_He"];
        C_ppt = props["C_ppt"];
        rho = props["rho"];
        a0 = props["a0"];
        d = props["d"];
        Z_i = props["Z_i"];
        r_ppt = props["r_ppt"];

        // int array_length = y.size();
        // dydt.resize(array_length);
        // const double &C_v = y[0];
        // const double &C_i = y[1];
        // const double &C_g = y[2];
        // const double &C_gv = y[3];
        // const double &C_2gv = y[4];
        // const double &C_2g = y[5];
        // const double &C_star = y[6];
        // const double &C_b = y[7];
        // const double &m = y[8];
        // const double &R = y[9];
        // const double &m_ppt = y[10];
        // const double &R_pt = y[11];
        // const double &M_gb = y[12];

        // Prepare time-dependent bubble-related parameters:
        std::tie(e3, epsilon, work) = get_bubble_props(y[9], y[8]);

        std::tie(e3_prime, epsilon_ppt, work_ppt) = get_bubble_props(y[11], y[10]);

        k_square = 4 * M_PI * y[9] * y[7] / Omega + rho;
        C_gb = a0 * a0 * sqrt(k_square) / (8 * d);

        // Sink concentrations
        bubble_sink = 4 * M_PI * y[9] * y[7] / Omega;
        precipitate_sink = 4 * M_PI * y[11] * C_ppt / Omega;

        C_s_v = (a0 * a0 / 48) * (rho + bubble_sink);
        C_s_i = (a0 * a0 / 48) * (Z_i * rho + bubble_sink);

        dydt[0] = (P + (beta * e1 + delta) * y[3] - (alpha * y[1] + beta * y[2] + gamma_val * (C_s_v + y[3] + 2 * (y[5] + y[4]) + 3 * y[6])) * y[0]);
        dydt[1] = P - alpha * (y[0] + y[3] + 2 * y[4] + 3 * y[6] + C_s_i) * y[1];
        dydt[2] = (G_He - beta * y[2] * (y[0] + 2 * y[2] + y[3] + 2 * y[5] + 2 * y[4] + C_gb + epsilon * y[7]) + delta * (y[3] + 2 * y[4] + 4 * y[5] + 3 * y[6] + y[8] * y[7] + y[12] + y[10] * C_ppt) + alpha * y[1] * y[3] + beta * (e1 * y[3] + e2 * y[4]));
        dydt[3] = (beta * y[2] * y[0] + (beta * e2 + 2 * delta) * y[4] - (beta * e1 + delta + alpha * y[1] + beta * y[2]) * y[3]);
        dydt[4] = (beta * y[2] * y[3] + 3 * delta * y[6] + 2 * gamma_val * y[0] * y[5] - (2 * beta * y[2] + 2 * delta + beta * e2 + 2 * alpha * y[1]) * y[4]);
        dydt[5] = (alpha * y[1] * y[4] + 2 * beta * y[2] * y[2] - (2 * delta + 2 * gamma_val * y[0] + 2 * beta * y[2]) * y[5]);
        dydt[6] = (2 * beta * y[2] * (y[4] + y[5]) - 3 * y[6] * (delta + alpha * y[1] + beta * y[2] + gamma_val * y[0]));
        dydt[7] = ((12 * beta * y[2] + 9 * gamma_val * y[0]) * y[6]) / y[8];
        dydt[8] = epsilon * beta * y[2] - delta * y[8];
        dydt[9] = (a0 * a0 / y[9]) * (gamma_val * y[0] - alpha * y[1] - gamma_val * (e3 - e4));
        if (dydt[9] < 0) dydt[9] = 0;
        dydt[10] = epsilon_ppt * beta * y[2] - delta * y[10];
        double r_p_equiv = sqrt(y[11] * y[11] + r_ppt * r_ppt);
        dydt[11] = (a0 * a0 / r_p_equiv) * (gamma_val * y[0] - alpha * y[1] - gamma_val * (e3_prime - e4));
        if (dydt[11] < 0) dydt[11] = 0;
        dydt[12] = beta * C_gb * y[2] - delta * y[12];

        // std::cout << "t: " << t << " dydt: ";
        // for (const auto& val : dydt) {
        //     std::cout << val << " ";
        // }
        // std::cout << std::endl;

        std::cout << "t: " << t << " y: ";
        for (const auto& val : y) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    };

template< class State , class Matrix >
void jacobian::operator()( const State &y , Matrix &J , const double &t , State &dfdt )
    {
        props = get_props();

        P = props["P"];
        alpha = props["alpha"];
        beta = props["beta"];
        gamma_val = props["gamma"];
        delta = props["delta"];
        e1 = props["e1"];
        e2 = props["e2"];
        e4 = props["e4"];
        e5 = props["e5"];
        Omega = props["Omega"];
        G_He = props["G_He"];
        C_ppt = props["C_ppt"];
        rho = props["rho"];
        a0 = props["a0"];
        d = props["d"];
        Z_i = props["Z_i"];
        r_ppt = props["r_ppt"];

        // int array_length = y.size();
        // dydt.resize(array_length);
        // const double &C_v = y[0];
        // const double &C_i = y[1];
        // const double &C_g = y[2];
        // const double &C_gv = y[3];
        // const double &C_2gv = y[4];
        // const double &C_2g = y[5];
        // const double &C_star = y[6];
        // const double &C_b = y[7];
        // const double &m = y[8];
        // const double &R = y[9];
        // const double &m_ppt = y[10];
        // const double &R_pt = y[11];
        // const double &M_gb = y[12];

        // Prepare time-dependent bubble-related parameters:
        std::tie(e3, epsilon, work) = get_bubble_props(y[9], y[8]);

        std::tie(e3_prime, epsilon_ppt, work_ppt) = get_bubble_props(y[11], y[10]);

        k_square = 4 * M_PI * y[9] * y[7] / Omega + rho;
        C_gb = a0 * a0 * sqrt(k_square) / (8 * d);

        // Sink concentrations
        bubble_sink = 4 * M_PI * y[9] * y[7] / Omega;
        precipitate_sink = 4 * M_PI * y[11] * C_ppt / Omega;

        C_s_v = (a0 * a0 / 48) * (rho + bubble_sink);
        C_s_i = (a0 * a0 / 48) * (Z_i * rho + bubble_sink);

        J(0, 0) = -(alpha * y[1] + beta * y[2] + gamma_val * (C_s_v + y[3] + 2 * (y[5] + y[4]) + 3 * y[6]));
        J(0, 1) = -(beta * y[0]);
        J(0, 2) = -(beta * y[0]);
        J(0, 3) = (beta * e1 + delta);
        J(0, 4) = 0;
        J(0, 5) = 0;
        J(0, 6) = 0;
        J(0, 7) = 0;
        J(0, 8) = 0;
        J(0, 9) = 0;
        J(0, 10) = 0;
        J(0, 11) = 0;
        J(0, 12) = 0;

        J(1, 0) = -(alpha * (y[0] + y[3] + 2 * y[4] + 3 * y[6] + C_s_i));
        J(1, 1) = -(alpha * (y[0] + y[3] + 2 * y[4] + 3 * y[6] + C_s_i));
        J(1, 2) = 0;
        J(1, 3) = 0;
        J(1, 4) = 0;
        J(1, 5) = 0;
        J(1, 6) = 0;
        J(1, 7) = 0;
        J(1, 8) = 0;
        J(1, 9) = 0;
        J(1, 10) = 0;
        J(1, 11) = 0;
        J(1, 12) = 0;

        J(2, 0) = -(beta * y[2] * (y[0] + 2 * y[2] + y[3] + 2 * y[5] + 2 * y[4] + C_gb + epsilon * y[7]));
        J(2, 1) = alpha * y[3];
        J(2, 2) = G_He - beta * (y[0] + 2 * y[2] + y[3] + 2 * y[5] + 2 * y[4] + C_gb + epsilon * y[7]);
        J(2, 3) = delta + alpha * y[1] + beta * y[2];
        J(2, 4) = beta * e1 + delta + alpha * y[1] + beta * y[2];
        J(2, 5) = 0;
        J(2, 6) = 0;
        J(2, 7) = beta * epsilon * y[2];
        J(2, 8) = 0;
        J(2, 9) = 0;
        J(2, 10) = delta * C_ppt;
        J(2, 11) = 0;
        J(2, 12) = beta * C_gb;

        J(3, 0) = beta * y[2];
        J(3, 1) = alpha * y[3];
        J(3, 2) = beta * y[0] + (beta * e2 + 2 * delta) * y[4];
        J(3, 3) = beta * y[0] + (beta * e2 + 2 * delta) * y[4];
        J(3, 4) = (beta * e1 + delta + alpha * y[1] + beta * y[2]);
        J(3, 5) = 0;
        J(3, 6) = 0;
        J(3, 7) = 0;
        J(3, 8) = 0;
        J(3, 9) = 0;
        J(3, 10) = 0;
        J(3, 11) = 0;
        J(3, 12) = 0;

        J(4, 0) = 2 * gamma_val * y[5];
        J(4, 1) = 0;
        J(4, 2) = beta * y[3] + 3 * delta * y[6];
        J(4, 3) = 0;
        J(4, 4) = beta * y[2] + 3 * delta * y[6];
        J(4, 5) = 2 * gamma_val * y[0] - (2 * beta * y[2] + 2 * delta + beta * e2 + 2 * alpha * y[1]);
        J(4, 6) = 3 * delta * y[4];
        J(4, 7) = 0;
        J(4, 8) = 0;
        J(4, 9) = 0;
        J(4, 10) = 0;
        J(4, 11) = 0;
        J(4, 12) = 0;

        J(5, 0) = 2 * beta * y[2];
        J(5, 1) = alpha * y[4];
        J(5, 2) = 4 * beta * y[2] + 2 * gamma_val * y[0] - (2 * delta + 2 * gamma_val * y[0] + 2 * beta * y[2]);
        J(5, 3) = 0;
        J(5, 4) = alpha * y[1] + 2 * beta * y[2] * y[2] - (2 * delta + 2 * gamma_val * y[0] + 2 * beta * y[2]);
        J(5, 5) = 0;
        J(5, 6) = 0;
        J(5, 7) = 0;
        J(5, 8) = 0;
        J(5, 9) = 0;
        J(5, 10) = 0;
        J(5, 11) = 0;
        J(5, 12) = 0;

        J(6, 0) = -3 * y[6] * gamma_val;
        J(6, 1) = -3 * y[6] * alpha;
        J(6, 2) = -3 * y[6] * beta;
        J(6, 3) = -3 * y[6] * gamma_val;
        J(6, 4) = -6 * y[6] * gamma_val;
        J(6, 5) = -6 * y[6] * gamma_val;
        J(6, 6) = -3 * (delta + alpha * y[1] + beta * y[2] + gamma_val * y[0]);
        J(6, 7) = 0;
        J(6, 8) = 0;
        J(6, 9) = 0;
        J(6, 10) = 0;
        J(6, 11) = 0;
        J(6, 12) = 0;

        J(7, 0) = 0;
        J(7, 1) = 0;
        J(7, 2) = 0;
        J(7, 3) = 0;
        J(7, 4) = 0;
        J(7, 5) = 0;
        J(7, 6) = 0;
        J(7, 7) = -delta;
        J(7, 8) = 0;
        J(7, 9) = 0;
        J(7, 10) = 0;
        J(7, 11) = 0;
        J(7, 12) = 0;

        J(8, 0) = 0;
        J(8, 1) = 0;
        J(8, 2) = 0;
        J(8, 3) = 0;
        J(8, 4) = 0;
        J(8, 5) = 0;
        J(8, 6) = 0;
        J(8, 7) = 0;
        J(8, 8) = -delta;
        J(8, 9) = 0;
        J(8, 10) = 0;
        J(8, 11) = 0;
        J(8, 12) = 0;

        J(9, 0) = 0;
        J(9, 1) = 0;
        J(9, 2) = 0;
        J(9, 3) = 0;
        J(9, 4) = 0;
        J(9, 5) = 0;
        J(9, 6) = 0;
        J(9, 7) = 0;
        J(9, 8) = 0;
        J(9, 9) = -(a0 * a0 / (y[9] * y[9]));
        J(9, 10) = 0;
        J(9, 11) = 0;
        J(9, 12) = 0;

        J(10, 0) = 0;
        J(10, 1) = 0;
        J(10, 2) = 0;
        J(10, 3) = 0;
        J(10, 4) = 0;
        J(10, 5) = 0;
        J(10, 6) = 0;
        J(10, 7) = 0;
        J(10, 8) = 0;
        J(10, 9) = 0;
        J(10, 10) = -delta;
        J(10, 11) = 0;
        J(10, 12) = 0;

        J(11, 0) = 0;
        J(11, 1) = 0;
        J(11, 2) = 0;
        J(11, 3) = 0;
        J(11, 4) = 0;
        J(11, 5) = 0;
        J(11, 6) = 0;
        J(11, 7) = 0;
        J(11, 8) = 0;
        J(11, 9) = 0;
        J(11, 10) = 0;
        J(11, 11) = -(a0 * a0 / (y[11] * y[11]));
        J(11, 12) = 0;

        J(12, 0) = 0;
        J(12, 1) = 0;
        J(12, 2) = 0;
        J(12, 3) = 0;
        J(12, 4) = 0;
        J(12, 5) = 0;
        J(12, 6) = 0;
        J(12, 7) = beta * C_gb;
        J(12, 8) = 0;
        J(12, 9) = 0;
        J(12, 10) = 0;
        J(12, 11) = 0;
        J(12, 12) = -delta;

        dfdt[0] = 0;
        dfdt[1] = 0;
        dfdt[2] = 0;
        dfdt[3] = 0;
        dfdt[4] = 0;
        dfdt[5] = 0;
        dfdt[6] = 0;
        dfdt[7] = 0;
        dfdt[8] = 0;
        dfdt[9] = 0;
        dfdt[10] = 0;
        dfdt[11] = 0;
        dfdt[12] = 0;
    };

// Explicit instantiation of the template for the specific types used
template void derivatives::operator()<boost::numeric::ublas::vector<double>>(const boost::numeric::ublas::vector<double> &y, boost::numeric::ublas::vector<double> &dydt, double t);
template void jacobian::operator()<boost::numeric::ublas::vector<double>, boost::numeric::ublas::matrix<double>>(const boost::numeric::ublas::vector<double> &y, boost::numeric::ublas::matrix<double> &J, const double &t, boost::numeric::ublas::vector<double> &dfdt);
// void derivatives(const state_type &y, state_type &dydt, double t) {
//     props = get_props();

//     double P = props["P"];
//     double alpha = props["alpha"];
//     double beta = props["beta"];
//     double gamma_val = props["gamma"];
//     double delta = props["delta"];
//     double e1 = props["e1"];
//     double e2 = props["e2"];
//     double e4 = props["e4"];
//     double e5 = props["e5"];
//     double Omega = props["Omega"];
//     double G_He = props["G_He"];
//     double C_ppt = props["C_ppt"];
//     double rho = props["rho"];
//     double a0 = props["a0"];
//     double d = props["d"];
//     double Z_i = props["Z_i"];
//     double r_ppt = props["r_ppt"];

//     int array_length = y.size();
//     dydt.resize(array_length);

//     const double &C_v = y[0];
//     const double &C_i = y[1];
//     const double &C_g = y[2];
//     const double &C_gv = y[3];
//     const double &C_2gv = y[4];
//     const double &C_2g = y[5];
//     const double &C_star = y[6];
//     const double &C_b = y[7];
//     const double &m = y[8];
//     const double &R = y[9];
//     const double &m_ppt = y[10];
//     const double &R_pt = y[11];
//     const double &M_gb = y[12];

//     // Prepare time-dependent bubble-related parameters:
//     double e3, epsilon, work;
//     std::tie(e3, epsilon, work) = get_bubble_props(R, m);

//     double e3_prime, epsilon_ppt, work_ppt;
//     std::tie(e3_prime, epsilon_ppt, work_ppt) = get_bubble_props(R_pt, m_ppt);

//     double k_square = 4 * M_PI * R * C_b / Omega + rho;
//     double C_gb = a0 * a0 * sqrt(k_square) / (8 * d);

//     // Sink concentrations
//     double bubble_sink = 4 * M_PI * R * C_b / Omega;
//     double precipitate_sink = 4 * M_PI * R_pt * C_ppt / Omega;

//     double C_s_v = (a0 * a0 / 48) * (rho + bubble_sink);
//     double C_s_i = (a0 * a0 / 48) * (Z_i * rho + bubble_sink);

//     // Rate equations
//     // dydt[0] = (P + (beta * e1 + delta) * C_gv - (alpha * C_i + beta * C_g + gamma * (C_s_v + C_gv + 2 * (C_2g + C_2gv) + 3 * C_star)) * C_v);
//     // dydt[1] = P - alpha * (C_v + C_gv + 2 * C_2gv + 3 * C_star + C_s_i) * C_i;
//     // dydt[2] = (G_He - beta * C_g * (C_v + 2 * C_g + C_gv + 2 * C_2g + 2 * C_2gv + C_gb + epsilon * C_b) + delta * (C_gv + 2 * C_2gv + 4 * C_2g + 3 * C_star + m * C_b + M_gb + m_ppt * C_ppt) + alpha * C_i * C_gv + beta * (e1 * C_gv + e2 * C_2gv));
//     // dydt[3] = (beta * C_g * C_v + (beta * e2 + 2 * delta) * C_2gv - (beta * e1 + delta + alpha * C_i + beta * C_g) * C_gv);
//     // dydt[4] = (beta * C_g * C_gv + 3 * delta * C_star + 2 * gamma * C_v * C_2g - (2 * beta * C_g + 2 * delta + beta * e2 + 2 * alpha * C_i) * C_2gv);
//     // dydt[5] = (alpha * C_i * C_2gv + 2 * beta * C_g * C_g - (2 * delta + 2 * gamma * C_v + 2 * beta * C_g) * C_2g);
//     // dydt[6] = (2 * beta * C_g * (C_2gv + C_2g) - 3 * C_star * (delta + alpha * C_i + beta * C_g + gamma * C_v));
//     // dydt[7] = ((12 * beta * C_g + 9 * gamma * C_v) * C_star) / m;
//     // dydt[8] = epsilon * beta * C_g - delta * m;
//     // dydt[9] = (a0 * a0 / R) * (gamma * C_v - alpha * C_i - gamma * (e3 - e4));
//     // if (dydt[9] < 0) dydt[9] = 0;
//     // dydt[10] = epsilon_ppt * beta * C_g - delta * m_ppt;
//     // double r_p_equiv = sqrt(R_pt * R_pt + r_ppt * r_ppt);
//     // dydt[11] = (a0 * a0 / r_p_equiv) * (gamma * C_v - alpha * C_i - gamma * (e3_prime - e4));
//     // if (dydt[11] < 0) dydt[11] = 0;
//     // dydt[12] = beta * C_gb * C_g - delta * M_gb;


//     dydt[0] = (P + (beta * e1 + delta) * y[3] - (alpha * y[1] + beta * y[2] + gamma * (C_s_v + y[3] + 2 * (y[5] + y[4]) + 3 * y[6])) * y[0]);
//     dydt[1] = P - alpha * (y[0] + y[3] + 2 * y[4] + 3 * y[6] + C_s_i) * y[1];
//     dydt[2] = (G_He - beta * y[2] * (y[0] + 2 * y[2] + y[3] + 2 * y[5] + 2 * y[4] + C_gb + epsilon * y[7]) + delta * (y[3] + 2 * y[4] + 4 * y[5] + 3 * y[6] + y[8] * y[7] + y[12] + y[10] * C_ppt) + alpha * y[1] * y[3] + beta * (e1 * y[3] + e2 * y[4]));
//     dydt[3] = (beta * y[2] * y[0] + (beta * e2 + 2 * delta) * y[4] - (beta * e1 + delta + alpha * y[1] + beta * y[2]) * y[3]);
//     dydt[4] = (beta * y[2] * y[3] + 3 * delta * y[6] + 2 * gamma * y[0] * y[5] - (2 * beta * y[2] + 2 * delta + beta * e2 + 2 * alpha * y[1]) * y[4]);
//     dydt[5] = (alpha * y[1] * y[4] + 2 * beta * y[2] * y[2] - (2 * delta + 2 * gamma * y[0] + 2 * beta * y[2]) * y[5]);
//     dydt[6] = (2 * beta * y[2] * (y[4] + y[5]) - 3 * y[6] * (delta + alpha * y[1] + beta * y[2] + gamma * y[0]));
//     dydt[7] = ((12 * beta * y[2] + 9 * gamma * y[0]) * y[6]) / y[8];
//     dydt[8] = epsilon * beta * y[2] - delta * y[8];
//     dydt[9] = (a0 * a0 / y[9]) * (gamma * y[0] - alpha * y[1] - gamma * (e3 - e4));
//     if (dydt[9] < 0) dydt[9] = 0;
//     dydt[10] = epsilon_ppt * beta * y[2] - delta * y[10];
//     double r_p_equiv = sqrt(y[11] * y[11] + r_ppt * r_ppt);
//     dydt[11] = (a0 * a0 / r_p_equiv) * (gamma * y[0] - alpha * y[1] - gamma * (e3_prime - e4));
//     if (dydt[11] < 0) dydt[11] = 0;
//     dydt[12] = beta * C_gb * y[2] - delta * y[12];
//     // Print intermediate values for debugging
//     std::cout.precision(10);
//     std::cout << "t: " << t << " dydt: ";
//     for (const auto& val : dydt) {
//         std::cout << val << " ";
//     }
//     std::cout << std::endl;
//     // std::cout << "dydt: ";
//     // for (const auto& val : dydt) {
//     //     std::cout << val << " ";
//     // }
//     // std::cout << std::endl;

// }

