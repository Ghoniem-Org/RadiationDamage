#include <iostream>
#include <boost/numeric/odeint.hpp>
#include "Utilities.h"
#include "Derivatives.h"

using namespace boost::numeric::odeint;
using namespace boost::numeric::ublas;

// Define the global variable
vector_type y;
int array_length;

int main() {
    // Read the Excel file
    read_excel("./material_data2.xlsx");

    // Initialize conditions
    // std::vector<double> y0;
    // double floor_val;
    // int array_length;
    // std::tie(y0, floor_val, array_length) = init_conditions();

    // initialize();
    // Define the initial state vector 'y'
    double floor_val = 1e-20;
    int array_length = 13;
    vector_type y(array_length,floor_val);
    y[8] = 2;
    y[9] = 5e-10;
    y[10] = 2;
    y[11] = 5e-10;

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
    
    // Solve the system of ODEs
    std::vector<vector_type> sol;
    std::vector<double> times;
    // Define the stepper
    // typedef runge_kutta4<std::vector<double>> stepper_type;
    // stepper_type stepper;

    // Integrate the system of ODEs with adaptive time steps based on t_eval
    // typedef runge_kutta_cash_karp54<state_type> error_stepper_type;
    // typedef controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
    // controlled_stepper_type controlled_stepper;

    // typedef adams_bashforth_moulton<state_type> stepper_type;
    // stepper_type stepper;
    // typedef bulirsch_stoer<state_type> stepper_type;
    // stepper_type stepper;

    // typedef implicit_midpoint<state_type> stepper_type;
    // stepper_type stepper;
    // typedef implicit_euler<state_type> stepper_type;
    // stepper_type stepper;


    typedef rosenbrock4<vector_type> stiff_stepper_type;
    stiff_stepper_type stepper;
    // typedef rosenbrock4_controller< stiff_stepper_type > controlled_stepper_type;
    // controlled_stepper_type stiff_stepper;

    // Integrate the system of ODEs
    size_t steps = integrate_times(make_controlled< rosenbrock4< double > >( 1.0e-20 , 1.5e-8 ) ,
            std::make_pair( derivatives() , jacobian() ) ,
                    y, t_eval.begin(), t_eval.end(), 1e-3, push_back_state_and_time(sol, times));

    // size_t steps = integrate_times(stepper,
    //                             derivatives,
    //                             y,
    //                             t_eval.begin(),
    //                             t_eval.end(),
    //                             1e-20,
    //                             push_back_state_and_time(sol, times)
    // );
    
    // Output the results
    for (size_t i = 0; i < sol.size(); ++i) {
        std::cout << "t: " << times[i] << " y: " << std::endl;
        for (const auto& val : sol[i]) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Total number of steps: " << steps << std::endl;

    return 0;
}
