# Project-3

To compile + link: 
g++ -larmadillo Particle.cpp PenningTrap.cpp Main.cpp -o2 -o Main.exe

To run: 
./Main.exe In dt tN

In is an integer between 1 and 7. dt and tN are floating point numbers specifying the timestep and the amount of microseconds that the simulation run for. 
In=1 simulates particle 1 with Forward Euler. In=2 simulates particle 2 with Forward Euler, In=3 simulates particle 1 with RK4, In = 4 simulates particles 1 and 2 
interacting with RK4, In=5 simulates particle 1 with Forward Euler for the given n's, In=6 simulates particle 1 with RK4 for the given n's, and In=7 simulates 20 
random particles. 

For In=1, 2, 3, 4 I have used dt=1e-4. For In = 5, 6, the timesteps are specified in the code so the dt input is irrelevant (but must stil be specified).

PenningTrap1.py has been run using "Run" in Spyder. 
