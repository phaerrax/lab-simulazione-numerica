#ifndef METROPOLIS_NORMAL
#define METROPOLIS_NORMAL

#include "metropolis.hh"

template <unsigned int dim>
class metropolis_normal : public metropolis<dim>
{
    public:
        metropolis_normal(double stdev);

        std::array<double, dim> step(const std::function<double (std::array<double, dim>)> &, Random &);

    private:
        double stdev;
};

#include "metropolis_normal.tcc"

#endif
