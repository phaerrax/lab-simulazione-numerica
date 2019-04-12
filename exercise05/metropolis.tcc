#include "metropolis.hh"

template <unsigned int dim>
metropolis<dim>::metropolis()
{
    proposals = 0;
    accepted_proposals = 0;

    return;
}

template <unsigned int dim>
double metropolis<dim>::get_acceptance_rate() const
{
    return static_cast<double>(accepted_proposals) / proposals;
}

template <unsigned int dim>
void metropolis<dim>::set_starting_point(const std::array<double, dim> & input_point)
{
    current_point = input_point;
    return;
}

template <unsigned int dim>
void metropolis<dim>::set_starting_point(std::array<double, dim> && input_point)
{
    current_point = input_point;
    return;
}
