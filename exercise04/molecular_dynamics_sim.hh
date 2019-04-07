#ifndef MOLECULAR_DYNAMICS_SIM
#define MOLECULAR_DYNAMICS_SIM

class molecular_dynamics_sim
{
    public:
        void input(const std::string & input_filename);
        // Collect from a file the initial configuration of the system.

        void initialise();
        // Generate an initial distribution of velocities, from which
        // simulate the position of the particles in the time step just
        // before the start.
        // The velocities are drawn from a Maxwell-Boltzmann distribution
        // in accord with the given temperature (this eliminates the need
        // to rescale the velocities after they are generated).
        // Knowing the position at the current and previous steps, and the
        // acceleration exerted on each particle, the algorithm can start.
        void initialise(const std::string &);
        // This "pre-initial" configuration can be also be supplied as a
        // text file; this way the position at the next step can be
        // calculated directly with the Verlet algorithm.
        // Then the velocity of the particles can be calculated, and it is
        // rescaled in order to match the previously set temperature.

        void move();
        // Advance the algorithm one step at a time with the Verlet method.

        // Read some thermodynamic properties of the system, calculated from
        // its current state.
        double potential_energy() const;
        double kinetic_energy() const;
        double temperature() const;
        double total_energy() const;

        void write_config(const std::string &) const;
        // Write the final configuration of the system in a file.
        void write_config_xyz(const std::string &, int) const;
        // Write the n-th configuration of the system in a file, in XYZ
        // format (to be read by ovito or similar programs).
        // This makes the visualisation of the progressive evolution of the
        // system as the algorithm advances possible.

    private:
        double force(unsigned int particle_index, unsigned int direction) const;
        // Calculate the force exerted on the particle at index
        // 'particle_index', in the direction 'direction'.

        double quotient(double) const;
        // Quotient the position/distance/etc. in order to reduce everything
        // to the unit cell of the lattice.
        // This imposes periodic boundary conditions on the unit cell of
        // the lattice.

        // Internal parameters
        double stima_pot, stima_kin, stima_etot, stima_temp;

        // averages
        double acc,att;

        //configuration
        const int m_part=108;
        double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
        double vx[m_part],vy[m_part],vz[m_part];

        // thermodynamical state
        int npart;
        double energy,temp,vol,rho,box,rcut;

        // simulation
        int nstep, iprint, seed;
        double delta;

};

#endif
