# include "Particle.hpp"

// Definitions of constructors
Particle::Particle(double q_in, double m_in, arma::vec r_in, arma::vec v_in)
{
  q = q_in; // Charge, e
  m = m_in; // Mass, u
  r = r_in; // Position, mu m
  v = v_in; // Velocity, (mu m)/(mu s)
}

