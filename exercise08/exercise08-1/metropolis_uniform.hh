#ifndef METROPOLIS_UNIFORM
#define METROPOLIS_UNIFORM

#include "metropolis.hh"

template <unsigned int dim>
class metropolis_uniform : public metropolis<dim>
{
    public:
        metropolis_uniform(double min, double max);

        std::array<double, dim> step(const std::function<double (std::array<double, dim>)> &, Random &);

    private:
        double min, max;
};

#include "metropolis_uniform.tcc"

#endif
