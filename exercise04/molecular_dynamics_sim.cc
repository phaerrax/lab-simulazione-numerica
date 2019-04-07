#include "molecular_dynamics_sim.hh"

void molecular_dynamics_sim::input(Random & rng)
{   
    // It is assumed that the random number generator is already initialised.

    std::string input_parameters_file("input.dat"),
                initial_config_file("config.0");

    double ep, ek, pr, et, vir; // ??

    std::cout << "Classic Lennard-Jones fluid\n
                  Molecular dynamics simulation in NVE ensemble\n\n
                  Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]\n\n
                  The program uses Lennard-Jones units" << std::endl;

    // Read the input parameters from a file.
    std::cout << "Read input parameters from file " << input_parameters_file << "." << std::endl;
    std::ifstream input_parameters(input_parameters_file);

    input_parameters >> temp;

    input_parameters >> npart;
    std::cout << "Number of particles = " << npart << std::endl;

    input_parameters >> rho;
    std::cout << "Density of particles = " << rho << std::endl;
    vol = static_cast<double>(npart)/rho;
    std::cout << "Volume of the simulation box = " << vol << std::endl;
    box = pow(vol,1.0/3.0);
    std::cout << "Edge of the simulation box = " << box << std::endl;

    input_parameters >> rcut;
    input_parameters >> delta;
    input_parameters >> nstep;
    input_parameters >> iprint;

    std::cout << "The program integrates Newton equations with the Verlet method, using\n
                  Time step = " << delta << ",\n
                  Number of steps = " << nstep << "." << std::endl;
    input_parameters.close(),

    // Read the initial configuration from a file.
    std::cout << "Read initial configuration from file " << initial_config_file << "." << std::endl;
    std::ifstream input_config(initial_config_file);

    for(unsigned int i = 0; i < npart; ++i)
    {
        input_config >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
    }
    input_config.close();

    // Simulate an initial velocity configuration.
    std::cout << "Prepare random velocities with center of mass velocity equal to zero " << std::endl << std::endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for(unsigned int i = 0; i < npart; ++i)
    {
        vx[i] = rng.Rannyu() - 0.5;
        vy[i] = rng.Rannyu() - 0.5;
        vz[i] = rng.Rannyu() - 0.5;

        sumv[0] += vx[i];
        sumv[1] += vy[i];
        sumv[2] += vz[i];
    }
    for(unsigned int idim = 0; idim < 3; ++idim)
        sumv[idim] /= static_cast<double>(npart);
    double sumv2 = 0.0, fs;
    for(unsigned int i = 0; i < npart; ++i)
    {
        vx[i] = vx[i] - sumv[0];
        vy[i] = vy[i] - sumv[1];
        vz[i] = vz[i] - sumv[2];
        sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= static_cast<double>(npart);

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    for(unsigned int i = 0; i < npart; ++i)
    {
        vx[i] *= fs;
        vy[i] *= fs;
        vz[i] *= fs;

        xold[i] = x[i] - vx[i] * delta;
        yold[i] = y[i] - vy[i] * delta;
        zold[i] = z[i] - vz[i] * delta;
    }
    return;
}

void molecular_dynamics_sim::move(void)
{ //Move particles with Verlet algorithm
    double xnew, ynew, znew, fx[npart], fy[npart], fz[npart];

    for(unsigned int i = 0; i < npart; ++i)
    { 
        // Calculate the force acting on particle i.
        fx[i] = force(i,0);
        fy[i] = force(i,1);
        fz[i] = force(i,2);
    }

    for(unsigned int i = 0; i < npart; ++i)
    { 
        // Verlet integration scheme
        xnew = quotient( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
        ynew = quotient( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
        znew = quotient( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

        vx[i] = quotient(xnew - xold[i])/(2.0 * delta);
        vy[i] = quotient(ynew - yold[i])/(2.0 * delta);
        vz[i] = quotient(znew - zold[i])/(2.0 * delta);

        xold[i] = x[i];
        yold[i] = y[i];
        zold[i] = z[i];

        x[i] = xnew;
        y[i] = ynew;
        z[i] = znew;
    }
    return;
}

double molecular_dynamics_sim::force(int ip, int idir)
{ //Compute forces as -Grad_ip V(r)
    double f=0.0;
    double dvec[3], dr;

    for (unsigned int i = 0; i < npart; ++i)
    {
        if(i != ip)
        {
            dvec[0] = quotient( x[ip] - x[i] );  // distance ip-i in pbc
            dvec[1] = quotient( y[ip] - y[i] );
            dvec[2] = quotient( z[ip] - z[i] );

            dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
            dr = sqrt(dr);
            if(dr < rcut)
            {
                f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
            }
        }
    }
    return f;
}

void molecular_dynamics_sim::measure()
{ //Properties measurement
    int bin;
    double v, t, vij;
    double dx, dy, dz, dr;
    std::ofstream Epot, Ekin, Etot, Temp;

    Epot.open("output_epot.dat",ios::app);
    Ekin.open("output_ekin.dat",ios::app);
    Temp.open("output_temp.dat",ios::app);
    Etot.open("output_etot.dat",ios::app);

    v = 0.0; //reset observables
    t = 0.0;

    //cycle over pairs of particles
    for(unsigned int i = 0; i < npart-1; ++i)
    {
        for(int j=i+1; j<npart; ++j)
        {
            dx = quotient( x[i] - x[j] );
            dy = quotient( y[i] - y[j] );
            dz = quotient( z[i] - z[j] );

            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);
            if(dr < rcut)
            {
                vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
                //Potential energy
                v += vij;
            }
        }          
    }

    //Kinetic energy
    for(unsigned int i = 0; i < npart; ++i)
        t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

    stima_pot = v/static_cast<double>(npart); //Potential energy
    stima_kin = t/static_cast<double>(npart); //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/static_cast<double>(npart); //Temperature
    stima_etot = (t+v)/static_cast<double>(npart); //Total enery

    Epot << stima_pot  << std::endl;
    Ekin << stima_kin  << std::endl;
    Temp << stima_temp << std::endl;
    Etot << stima_etot << std::endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}

void molecular_dynamics_sim::write_config(const std::string & output_file) const
{ 
    // Write the final configuration on file.
    std::cout << "Print final configuration to file " << output_file << "." << std::endl;
    std::ofstream output_file(output_file);

    for(unsigned int i = 0; i < npart; ++i)
    {
        output_file << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << std::endl;
    }
    output_file.close();
    return;
}

void molecular_dynamics_sim::write_config_xyz(const std::string & output_file, int nconf) const
{ 
    // Write configuration in .xyz format.
    std::ofstream output_file;

    output_file.open("frames/config_" + std::to_string(nconf) + ".xyz");
    output_file << npart << std::endl;
    output_file << "This is only a comment!" << std::endl;
    for(unsigned int i = 0; i < npart; ++i)
    {
        output_file << "LJ  " << quotient(x[i]) << "   " <<  quotient(y[i]) << "   " << quotient(z[i]) << std::endl;
    }
    output_file.close();
}

double molecular_dynamics_sim::quotient(double r) const
{  
    return r - box * std::round(r / box);
}
