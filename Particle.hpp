#ifndef __Particle_hpp__
#define __Particle_hpp__

# include <iostream>
# include <armadillo>

class Particle 
{
    public:
        double q;
        double m;
        arma::vec r;
        arma::vec v;
        Particle(double q_in, double m_in, arma::vec r_in, arma::vec v_in);

};


#endif