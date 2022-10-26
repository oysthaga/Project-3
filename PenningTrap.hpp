#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

# include <iostream>
# include <vector>
# include <armadillo>
# include <cmath>

class PenningTrap 
{
    public:
      double B0;
      double V0;
      double d;
      double V0dd;
      double k_e;
      std::vector<Particle> Particles;

      // Constructor
      PenningTrap(double B0_in, double V0_in, double d_in, double V0dd_in, double k_e_in, std::vector<Particle> Particles_in);

      // Add a particle to the trap
      void add_particle(Particle p_in);

  // External electric field at point r=(x,y,z)
  arma::vec external_E_field(arma::vec r, double t=0, double omV=0, double f=0);  

  // External magnetic field at point r=(x,y,z)
  arma::vec external_B_field(arma::vec r);  

  // Force on particle_i from particle_j
  arma::vec force_particle(int i, int j);

  // The total force on particle_i from the external fields
  arma::vec total_force_external(int i, double t=0, double omV=0, double f=0);

  // The total force on particle_i from the other particles
  arma::vec total_force_particles(int i);

  // The total force on particle_i from both external fields and other particles
  arma::vec total_force(int i, double t=0, bool CoulombOn=true, double omV=0, double f=0);

  // Evolve the system one time step (dt) using Runge-Kutta 4th order
  void evolve_RK4(double dt, double t0=0, bool CoulombOn=true, double omV=0, double f=0);

  // Evolve the system one time step (dt) using Forward Euler
  void evolve_forward_Euler(double dt, double t=0, bool CoulombOn=true, double omV=0, double f=0);

  // Count particles remaining in trap 
  int Particle_counter();
};


#endif