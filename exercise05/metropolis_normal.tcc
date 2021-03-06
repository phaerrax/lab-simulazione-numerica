#include "metropolis_normal.hh"
#include <stdexcept>

template <unsigned int dim>
metropolis_normal<dim>::metropolis_normal(double istdev):
    metropolis<dim>()
{
    if(istdev <= 0)
		throw std::invalid_argument("Standard deviation must be positive");
	stdev = istdev;
    return;
}

template <unsigned int dim>
std::array<double, dim> metropolis_normal<dim>::step(
        const std::function<double (std::array<double, dim>)> & f,
        Random & rng
        )
{
    // Base classes members are hidden to this class, since metropolis
    // and this one are templated.
    // In order to use those members, they have to be "manually" brought 
    // into scope, or used with the 'this' pointer.

    // Generate an uniformly generated new point.
    std::array<double, dim> proposed_point;
    for(unsigned int i = 0; i < dim; ++i)
        proposed_point[i] = this->current_point[i] + rng.Gauss(0, stdev);
    ++this->proposals;

    double ratio = f(proposed_point) / f(this->current_point);
    // Accept the new point with probability min(1, ratio): if ratio is
    // already greater than one we can avoid generating a random number.
    if(ratio >= 1)
    {
        this->current_point = proposed_point;
        ++this->accepted_proposals;
    }
    else
        // Draw a random number and accept the proposal if that number
        // is less than ratio.
        if(rng.Rannyu() < ratio)
        {
            this->current_point = proposed_point;
            ++this->accepted_proposals;
        }

    return this->current_point;
}
