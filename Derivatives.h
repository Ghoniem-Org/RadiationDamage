#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include <vector>
#include "Utilities.h"

// Declare functions
// void derivatives(const state_type &y, state_type &dydt, double t);
extern vector_type y;
extern int array_length;

template< class State >
void load_variables(const State &y);

struct derivatives
{
    template< class State >
    void operator()(const State &y , State &dydt , double t );
};

struct jacobian
{
    template< class State , class Matrix >
    void operator()( const State &y , Matrix &J , const double &t , State &dfdt );
};

// The observer function
struct push_back_state_and_time
{
    std::vector<vector_type> &m_states;
    std::vector<double> &m_times;

    push_back_state_and_time(std::vector<vector_type> &states, std::vector<double> &times)
            : m_states(states), m_times(times) {}

    void operator()(const vector_type &y, double t)
    {
        m_states.push_back(y);
        m_times.push_back(t);
    }
};
// struct push_back_state_and_time
// {
//     std::vector<state_type> &m_states;
//     std::vector<double> &m_times;

//     push_back_state_and_time(std::vector<state_type> &states, std::vector<double> &times)
//             : m_states(states), m_times(times) {}

//     void operator()(const state_type &y, double t)
//     {
//         m_states.push_back(y);
//         m_times.push_back(t);
//     }
// };
#endif // DERIVATIVES_H