#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <map>
#include <vector>
#include <tuple>


#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/odeint.hpp>

// using namespace boost::numeric::odeint;
// using namespace boost::numeric::ublas;

// Global variables
// Define state type
// typedef std::vector<double> state_type;
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;
// typedef boost::numeric::ublas::vector<double> state_type;
extern vector_type y;
extern double floor_val;
extern int array_length;

// Declare functions
// std::tuple<std::vector<double>, double, int> init_conditions();
void initialize();
void read_excel(const std::string& file_path);
std::map<std::string, double> get_props();
std::tuple<double, double, double> get_bubble_props(double R, double m);

#endif // UTILITIES_H
