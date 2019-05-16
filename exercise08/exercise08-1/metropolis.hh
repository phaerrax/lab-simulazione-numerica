#ifndef METROPOLIS
#define METROPOLIS

#include <array>
#include <functional>
#include "random.hh"

template <unsigned int dim>
class metropolis
{
    public:
        metropolis();

        virtual std::array<double, dim> step(const std::function<double (std::array<double, dim>)> &, Random &) = 0;
        // Propose a new step in the chain and accept it with a given
        // probability. If it is a success, returns the new step, otherwise
        // returns the current step.

        void set_starting_point(const std::array<double, dim> &);
        void set_starting_point(std::array<double, dim> &&);
        // Set the current point to the given point.

        double get_acceptance_rate() const;
        // For debug purposes, returns the ratio of accepted steps vs the
        // total number of proposed steps.

    protected:
        unsigned int accepted_proposals,
                     proposals;

        std::array<double, dim> current_point;
};

#include "metropolis.tcc"

#endif
