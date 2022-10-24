# include "Particle.hpp"
# include "PenningTrap.hpp"


// Definitions of constructors
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, double V0dd_in, double k_e_in, std::vector<Particle> Particles_in, bool Vt_in)
{

    B0 = B0_in;
    V0 = V0_in; 
    d = d_in; 
    V0dd = V0dd_in;
    k_e = k_e_in;
    Particles = Particles_in;
    Vt = Vt_in;
}

// Add a particle to the trap
void PenningTrap::add_particle(Particle p_in)
{
    Particles.push_back(p_in);
}



// External electric field at point r=(x,y,z)
arma::vec PenningTrap::external_E_field(arma::vec r, double t, double omV, double f)
{
    if ( sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]) > d )
    { 
        return arma::vec{0, 0, 0};
    }
    if ( Vt==true )
    {
        return (-V0dd*(1 + f*cos(omV * t))*arma::vec{ -r(0), -r(1), 2*r(2)} ); 
    }
    return (-V0dd*arma::vec{ -r(0), -r(1), 2*r(2)} );
} 

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap::external_B_field(arma::vec r)
{
    if ( sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]) > d )
    { 
        return arma::vec{0, 0, 0};
    }
    return ( arma::vec{0, 0, B0} );
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j)
{
    arma::vec r_ij = Particles[i].r-Particles[j].r;
    double norm_r_ij = sqrt(r_ij[0]*r_ij[0]+r_ij[1]*r_ij[1]+r_ij[2]*r_ij[2]);
    
    //return arma::vec{0., 0., 0.};

    double norm_r_ij_3 =pow(norm_r_ij, 3);
    return ( k_e*Particles[i].q*Particles[j].q*(r_ij)/norm_r_ij_3 );
}


// The total force on particle_i from the external fields
arma::vec PenningTrap::total_force_external(int i, double t, double omV, double f)
{
    return ( Particles[i].q*external_E_field(Particles[i].r, t, omV, f) + Particles[i].q*( arma::cross(Particles[i].v, external_B_field(Particles[i].r)) ) );
}

// The total force on particle_i from the other particles
arma::vec PenningTrap::total_force_particles(int i)
{
    arma::vec F = arma::vec{0,0,0};
    
    for (int j = 0; j < Particles.size(); j++)
    {
        if (j != i)
        {
            F += force_particle(i, j);
        }
    }
    return (F);
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap::total_force(int i, double t, bool CoulombOn, double omV, double f)
{
    if (CoulombOn=true)
    {
        return ( total_force_external(i, t) + total_force_particles(i) );
    }
    else
    {
        return ( total_force_external(i, t) ); 
    }
}

// Evolve the system one time step (dt) using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt, double t0, bool CoulombOn, double omV, double f)
{
    std::vector<arma::vec> K1r;
    std::vector<arma::vec> K1v;
    std::vector<arma::vec> K2r;
    std::vector<arma::vec> K2v;
    std::vector<arma::vec> K3r;
    std::vector<arma::vec> K3v;
    std::vector<arma::vec> K4r;
    std::vector<arma::vec> K4v;

    std::vector<Particle> p = Particles;

    for (int j = 0; j < p.size(); j++)
    {
        K1r.push_back( p[j].v*dt );
        K1v.push_back( (total_force(j, t0, CoulombOn, omV, f)/p[j].m)*dt );
        p[j].r = Particles[j].r + K1r[j]/2;
        p[j].v = Particles[j].v +  K1v[j]/2;
    }
    for (int j = 0; j < p.size(); j++)
    {
        K2r.push_back( p[j].v*dt );
        K2v.push_back( (total_force(j, t0+dt/2, CoulombOn, omV, f)/p[j].m)*dt );
        p[j].r = Particles[j].r + K2r[j]/2;
        p[j].v = Particles[j].v + K2v[j]/2;
    }
    for (int j = 0; j < p.size(); j++)
    {
        K3r.push_back( p[j].v*dt );
        K3v.push_back( (total_force(j, t0+dt, CoulombOn, omV, f)/p[j].m)*dt );
        p[j].r = Particles[j].r + K3r[j];
        p[j].v = Particles[j].v + K3r[j];
    }
    for (int j = 0; j < p.size(); j++)
    {
        K4r.push_back( p[j].v*dt );
        K4v.push_back( (total_force(j, t0+dt, CoulombOn, omV, f)/p[j].m)*dt );
        Particles[j].r = Particles[j].r + (K1r[j] + 2*K2r[j] + 2*K3r[j] + K4r[j])/6;
        Particles[j].v = Particles[j].v + (K1v[j] + 2*K2v[j] + 2*K3v[j] + K4v[j])/6;
    }
}

// Evolve the system one time step (dt) using Forward Euler
void PenningTrap::evolve_forward_Euler(double dt, double t0, bool CoulombOn, double omV, double f)
{
    for (int j = 0; j < Particles.size(); j++)
    {
        Particles[j].v = Particles[j].v + (total_force(j, t0)/Particles[j].m, CoulombOn, omV, f)*dt;
        Particles[j].r = Particles[j].r + Particles[j].v*dt;
    }
}

int PenningTrap::Particle_counter()
{
    int count = 0;
    for (int i = 0; i < Particles.size(); i++)
    {
        double norm_r = sqrt(Particles[i].r[0]*Particles[i].r[0]+Particles[i].r[1]*Particles[i].r[1]+Particles[i].r[2]*Particles[i].r[2]);
        if (norm_r < d)
        {
            count += 1;
        }
        
    }
    return count;
}
