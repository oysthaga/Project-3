# include "Particle.hpp"
# include "PenningTrap.hpp"
//# include <string>
int main(int argc, char* argv[])
{

    // Check number of command-line arguments
    if (argc != 4)  // Expects 4 command-line arguments
    {
        // Get the name of the executable file
        std::string executable_name = argv[0];
        std::cerr << "Error: Wrong number of input arguments." << std::endl;
        std::cerr << "Usage: " << executable_name << " <some integer> " << std::endl;
        // Exit program with non-zero return code to indicate a problem
        return 1;   
    }

    int in = atoi(argv[1]);
    double dt = atof(argv[2]);
    double tN =atof(argv[3]);
    double t0 = 0.;

    double B0 = 9.65e+1; // u/( (mu s)e )
    double V0 = 2.41e+6; // u(mu m)^2 / ( (mu s)^2 e) 
    double d = 500; // mu m 
    double k_e = 1.38935333e+5; // u(mu m)^3/( (mu s)^2 e^2 )

    double q = 1.; // e
    double m = 40.078; // u

    arma::vec r10 = arma::vec{20.,  0., 20.}; // mu m
    arma::vec v10 = arma::vec{ 0., 25.,  0.}; // (mu m)/(mu s) 
    arma::vec r20 = arma::vec{25., 25.,  0.}; // mu m
    arma::vec v20 = arma::vec{ 0., 40.,  5.}; // (mu m)/(mu s)
    Particle P1 = Particle(q, m, r10, v10);
    Particle P2 = Particle(q, m, r20, v20);

    std::vector<Particle> Particles;
    PenningTrap obj = PenningTrap(B0, V0, d, V0/(d*d), k_e, Particles);


    if (in == 1) // Simulate P1 with Forward Euler
    {
        obj.add_particle(P1);
        arma::vec time = arma::regspace(t0, dt, tN);
        int N = time.size();
        arma::mat v1_Eu = arma::mat(3, N);
        arma::mat r1_Eu = arma::mat(3, N);

        v1_Eu.col(0) = v10;
        r1_Eu.col(0) = r10;
        for(int i=1; i < N; i++) 
        {
            obj.evolve_forward_Euler(dt);
            v1_Eu.col(i) = obj.Particles[0].v; 
            r1_Eu.col(i) = obj.Particles[0].r; 
        }

        v1_Eu.save("v1Eu.bin");
        r1_Eu.save("r1Eu.bin"); 
    }

    if (in == 2) // Simulate P1 with RK4
    {
        obj.add_particle(P1);
        arma::vec time = arma::regspace(t0, dt, tN);
        int N = time.size();
        arma::mat v1_RK = arma::mat(3, N);
        arma::mat r1_RK = arma::mat(3, N);

        v1_RK.col(0) = v10;
        r1_RK.col(0) = r10;
        for(int i=1; i < N; i++)
        {
            obj.evolve_RK4(dt);
            v1_RK.col(i) = obj.Particles[0].v; 
            r1_RK.col(i) = obj.Particles[0].r; 
        }
    
    
        v1_RK.save("v1RK.bin");
        r1_RK.save("r1RK.bin"); 
    }

    if (in == 3) // Simulate P2 with RK4
    {
        obj.add_particle(P2);
        arma::vec time = arma::regspace(t0, dt, tN);
        int N = time.size();
        arma::mat v2_RK = arma::mat(3, N);
        arma::mat r2_RK = arma::mat(3, N);

        v2_RK.col(0) = v20;
        r2_RK.col(0) = r20;
        for(int i=1; i < N; i++)
        {
            obj.evolve_RK4(dt);
            v2_RK.col(i) = obj.Particles[0].v; 
            r2_RK.col(i) = obj.Particles[0].r; 
        }
    
    
        v2_RK.save("v2RK.bin");
        r2_RK.save("r2RK.bin"); 
    }

    if (in == 4) // Simulate P1 and P2 (interacting) with RK4
    {
        obj.add_particle(P1);
        obj.add_particle(P2);
        arma::vec time = arma::regspace(t0, dt, tN);
        int N = time.size();
        arma::mat v1_RK = arma::mat(3, N);
        arma::mat r1_RK = arma::mat(3, N);
        arma::mat v2_RK = arma::mat(3, N);
        arma::mat r2_RK = arma::mat(3, N);

        v1_RK.col(0) = v10;
        r1_RK.col(0) = r10;
        v2_RK.col(0) = v20;
        r2_RK.col(0) = r20;

        for(int i=1; i < N; i++)
        {
            obj.evolve_RK4(dt);
            v1_RK.col(i) = obj.Particles[0].v; 
            r1_RK.col(i) = obj.Particles[0].r; 
            v2_RK.col(i) = obj.Particles[1].v; 
            r2_RK.col(i) = obj.Particles[1].r; 
        }
    
    
        v1_RK.save("v1RK_interacting.bin");
        r1_RK.save("r1RK_interacting.bin"); 
        v2_RK.save("v2RK_interacting.bin");
        r2_RK.save("r2RK_interacting.bin"); 
    }

    if (in == 5) // Simulate P1 with Forward Euler
    {
        std::vector<int> n = {4000, 8000, 16000, 32000};
        for (int i = 0; i < n.size(); i++) 
        {
            double dt_n = 50./n[i]; // mu s 
            std::vector<Particle> Particles_n;
            PenningTrap obj_n = PenningTrap(B0, V0, d, V0/(d*d), k_e, Particles_n);
            obj_n.add_particle(P1);
            arma::vec time = arma::regspace(t0, dt_n, tN);
            int N = time.size();
            arma::mat r1_Eu = arma::mat(3, N);

            r1_Eu.col(0) = r10;
            for(int i=1; i < N; i++) 
            {
                obj_n.evolve_forward_Euler(dt_n);
                r1_Eu.col(i) = obj_n.Particles[0].r; 
            }
            std::string rname = "r1_Eu_";
            rname += std::to_string(i);
            rname += ".bin";

            r1_Eu.save(rname); 
        }
    }

    if (in == 6) // Simulate P1 with RK4
    {
        std::vector<int> n = {4000, 8000, 16000, 32000};
        for (int i = 0; i < n.size(); i++) 
        {
            double dt_n = 50./n[i]; // mu s 
            std::vector<Particle> Particles_n;
            PenningTrap obj_n = PenningTrap(B0, V0, d, V0/(d*d), k_e, Particles_n);
            obj_n.add_particle(P1);
            arma::vec time = arma::regspace(t0, dt_n, tN);
            int N = time.size();
            arma::mat r1_RK = arma::mat(3, N);

            r1_RK.col(0) = r10;
            for(int i=1; i < N; i++) 
            {
                obj_n.evolve_RK4(dt_n);
                r1_RK.col(i) = obj_n.Particles[0].r; 
            }
            std::string rname = "r1_RK_";
            rname += std::to_string(i);
            rname += ".bin";

            r1_RK.save(rname); 
        }
    }

    if (in == 7)
    {
        // Construct a Mersenne Twister 19937 random number generator with a given seed
        std::mt19937 generator(1415184231);
        // Construct a distribution object for the uniform distribution on [0,1)
        std::normal_distribution<double> distribution(0.0 ,1.0);

        std::vector<double> f = {0.1, 0.4, 0.7};
        for  (int i = 0; i < f.size(); i++)
        {
            PenningTrap obj = PenningTrap(B0, V0, d, V0/(d*d), k_e, Particles, true);
            
            arma::vec time = arma::regspace(t0, dt, tN);
            arma::vec omV  = arma::regspace(0.2, 0.02, 2.5); // MHz
            int N =  omV.size();
            int M = time.size();
/*
            std::cout << "f =" << f[i];
            std::cout <<"\n";
*/
            for (int j = 0; j<20; j++)
            {


                arma::vec rj = arma::vec{distribution(generator), distribution(generator), distribution(generator)} * 0.1 * obj.d;
                arma::vec vj = arma::vec{distribution(generator), distribution(generator), distribution(generator)} * 0.1 * obj.d;
                obj.add_particle( Particle(q, m, rj, vj) ); // Add a particle with random initial velocity and position
/*
                std::cout << rj;
                std::cout <<"\n";
                std::cout << vj;
                std::cout <<"\n";
*/
            }
           
            arma::vec NumParticles = arma::vec(N);
            for(int j=0; j < N; j++) 
                {
                    for(int k=1; k < M; k++)  
                    {
                        obj.evolve_RK4(dt, time[k], false, omV[j], f[i]);
                    }
                    NumParticles(j) = obj.Particle_counter();
                }
            std::string name = "NumParticles";
            name += std::to_string(i);
            name += ".bin";

            NumParticles.save(name); 
         
        }
    }

    return 0;
}