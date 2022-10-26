import pyarma as pa
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

"""
v1Eu = pa.mat()
r1Eu = pa.mat()
v1Eu.load("v1Eu.bin")
r1Eu.load("r1Eu.bin")

v1RK = pa.mat()
r1RK = pa.mat()
v1RK.load("v1RK.bin")
r1RK.load("r1RK.bin")


v2RK = pa.mat()
r2RK = pa.mat()
v2RK.load("v2RK.bin")
r2RK.load("r2RK.bin")


v1RK_interacting = pa.mat()
r1RK_interacting = pa.mat()
v1RK_interacting.load("v1RK_interacting.bin")
r1RK_interacting.load("r1RK_interacting.bin")

v2RK_interacting = pa.mat()
r2RK_interacting = pa.mat()
v2RK_interacting.load("v2RK_interacting.bin")
r2RK_interacting.load("r2RK_interacting.bin")

"""

r1Eu0 = pa.mat()
r1Eu1 = pa.mat()
r1Eu2 = pa.mat()
r1Eu3 = pa.mat()
r1Eu0.load("r1_Eu_0.bin")
r1Eu1.load("r1_Eu_1.bin")
r1Eu2.load("r1_Eu_2.bin")
r1Eu3.load("r1_Eu_3.bin")
r1RK0 = pa.mat()
r1RK1 = pa.mat()
r1RK2 = pa.mat()
r1RK3 = pa.mat()
r1RK0.load("r1_RK_0.bin")
r1RK1.load("r1_RK_1.bin")
r1RK2.load("r1_RK_2.bin")
r1RK3.load("r1_RK_3.bin")


# COMPARING THE ANALYTICAL AND NUMERICAL SOLUTIONS 
q = 1 # e
V0 = 2.41e+6 # u(mu m)^2 / ( (mu s)^2 e) 
B0 = 9.65e+1 # u/( (mu s)e )
m = 40.078 # u
d = 500 # mu m 

x0 = 20. # mu m
znull = 20. # mu m
v0 = 25. # mu m / mu s
om0 = (q*B0)/m
omz = np.sqrt( (2*q*V0)/(m*d**2) )
omplus = ( om0 + np.sqrt(om0**2-2*omz**2) )/2
omminus =  ( om0 - np.sqrt(om0**2-2*omz**2) )/2
phiplus = 0
phiminus = 0
Aplus = (v0 + omminus*x0)/(omminus - omplus)
Aminus = -(v0 + omplus*x0)/(omminus - omplus)


t0 = np.linspace(0,50,np.shape(r1Eu0)[1]) 
t1 = np.linspace(0,50,np.shape(r1Eu1)[1])
t2 = np.linspace(0,50,np.shape(r1Eu2)[1])
t3 = np.linspace(0,50,np.shape(r1Eu3)[1])
z0 = znull*np.cos(omz*t0)
z1 = znull*np.cos(omz*t1)
z2 = znull*np.cos(omz*t2)
z3 = znull*np.cos(omz*t3)

f0 = Aplus*np.exp(-1j*(omplus*t0+phiplus)) + Aminus*np.exp(-1j*(omminus*t0+phiminus))
f1 = Aplus*np.exp(-1j*(omplus*t1+phiplus)) + Aminus*np.exp(-1j*(omminus*t1+phiminus))
f2 = Aplus*np.exp(-1j*(omplus*t2+phiplus)) + Aminus*np.exp(-1j*(omminus*t2+phiminus))
f3 = Aplus*np.exp(-1j*(omplus*t3+phiplus)) + Aminus*np.exp(-1j*(omminus*t3+phiminus)) 

r0 = np.array([f0.real, f0.imag, z0]); r1 = np.array([f1.real, f1.imag, z1])
r2 = np.array([f2.real, f2.imag, z2]); r3 = np.array([f3.real, f3.imag, z3])

"""
errEu0 = np.linalg.norm((r0 - r1Eu0)/r0, axis=0)
errEu1 = np.linalg.norm((r1 - r1Eu1)/r1, axis=0) 
errEu2 = np.linalg.norm((r2 - r1Eu2)/r2, axis=0)
errEu3 = np.linalg.norm((r3 - r1Eu3)/r3, axis=0)
errRK0 = np.linalg.norm((r0 - r1RK0)/r0, axis=0)
errRK1 = np.linalg.norm((r1 - r1RK1)/r1, axis=0)
errRK2 = np.linalg.norm((r2 - r1RK2)/r2, axis=0)
errRK3 = np.linalg.norm((r3 - r1RK3)/r3, axis=0)
"""
errEu0 = np.linalg.norm(r0 - r1Eu0, axis=0) / np.linalg.norm(r0, axis=0)
errEu1 = np.linalg.norm(r1 - r1Eu1, axis=0) / np.linalg.norm(r1, axis=0)
errEu2 = np.linalg.norm(r2 - r1Eu2, axis=0) / np.linalg.norm(r2, axis=0)
errEu3 = np.linalg.norm(r3 - r1Eu3, axis=0) / np.linalg.norm(r3, axis=0)
errRK0 = np.linalg.norm(r0 - r1Eu0, axis=0) / np.linalg.norm(r0, axis=0)
errRK1 = np.linalg.norm(r1 - r1Eu1, axis=0) / np.linalg.norm(r1, axis=0)
errRK2 = np.linalg.norm(r2 - r1Eu2, axis=0) / np.linalg.norm(r2, axis=0)
errRK3 = np.linalg.norm(r3 - r1Eu3, axis=0) / np.linalg.norm(r3, axis=0)

DeltaEu0 = max(np.linalg.norm(r0 - r1Eu0, axis=0)) 
DeltaEu1 = max(np.linalg.norm(r1 - r1Eu1, axis=0))
DeltaEu2 = max(np.linalg.norm(r2 - r1Eu2, axis=0)) 
DeltaEu3 = max(np.linalg.norm(r3 - r1Eu3, axis=0))
DeltaRK0 = max(np.linalg.norm(r0 - r1RK0, axis=0)) 
DeltaRK1 = max(np.linalg.norm(r1 - r1RK1, axis=0))
DeltaRK2 = max(np.linalg.norm(r2 - r1RK2, axis=0)) 
DeltaRK3 = max(np.linalg.norm(r3 - r1RK3, axis=0))

# h_k/h_(k-1) = (n_(k-1))/n_k = 2 
rerrEu = (1/3)*(np.log(DeltaEu1/DeltaEu0)+ np.log(DeltaEu2/DeltaEu1) 
                + np.log(DeltaEu3/DeltaEu2))/np.log(2) 
rerrRK = (1/3)*(np.log(DeltaRK1/DeltaRK0) + np.log(DeltaRK2/DeltaRK1) 
                + np.log(DeltaRK3/DeltaRK2))/np.log(2) 

print(f'r_(err, Euler) = {rerrEu}')
print(f'r_(err, RK4) = {rerrRK}')

"""
# PARTICLES LEFT
NumParticles0 = pa.mat()
NumParticles1 = pa.mat()
NumParticles2 = pa.mat()

NumParticles0.load("NumParticles0.bin")
NumParticles1.load("NumParticles1.bin")
NumParticles2.load("NumParticles2.bin")


omV = np.linspace(0.2,2.5,len(NumParticles0))

# PLOT z AS A FUNCTION OF TIME
N1 = np.shape(r1RK)[1]
t = np.linspace(0,50,N1)
plt.figure()
plt.plot(t, r1Eu[2,:], label='Euler')
plt.plot(t, r1RK[2,:], label='RK4')
plt.xlabel('$t [\mu s]$'); plt.ylabel('$z [\mu m]$')
plt.legend()
plt.show()

# PLOT xy-TRAJECTORY WITHOUT INTERACTION
plt.figure()
plt.title("Without interaction")
plt.plot(r1RK[0,:], r1RK[1,:], label = "$P1$")
plt.plot(r2RK[0,:], r2RK[1,:], label = "$P2$")
plt.xlabel("$x [\mu m]$"); plt.ylabel("$y [\mu m]$")
plt.legend()
plt.show()


# PLOT xy-TRAJECTORY WITH INTERACTION
plt.figure()
plt.title("With interaction")
plt.plot(r1RK_interacting[0,:], r1RK_interacting[1,:], label = "$P1$")
plt.plot(r2RK_interacting[0,:], r2RK_interacting[1,:], label = "$P2$")
plt.xlabel("$x [\mu m]$"); plt.ylabel("$y [\mu m]$")
plt.legend()
plt.axis('equal')
plt.show()


# PLOT v_x AS A FUNCTION OF x WITHOUT INTERACTION
plt.figure()
plt.title("Without interaction")
plt.plot(r1RK[0,:], v1RK[0,:], label = "$P1$")
plt.plot(r2RK[0,:], v2RK[0,:], label = "$P2$")
plt.xlabel("$x [\mu m]$"); plt.ylabel("$v_x [\mu m/ \mu s]$")
plt.legend()
plt.show()

# PLOT v_z AS A FUNCTION OF z WITHOUT INTERACTION
plt.figure()
plt.title("Without interaction")
plt.plot(r1RK[2,:], v1RK[2,:], label = "$P1$")
plt.plot(r2RK[2,:], v2RK[2,:], label = "$P2$")
plt.xlabel("$z [\mu m]$"); plt.ylabel("$v_z [\mu m/ \mu s]$")
plt.legend()
plt.show()

# PLOT v_x AS A FUNCTION OF x WITH INTERACTION
plt.figure()
plt.title("With interaction")
plt.plot(r1RK_interacting[0,:], v1RK_interacting[0,:], label = "$P1$")
plt.plot(r2RK_interacting[0,:], v2RK_interacting[0,:], label = "$P2$")
plt.xlabel("$x [\mu m]$"); plt.ylabel("$v_x [\mu m/ \mu s]$")
plt.legend()
plt.show()

# PLOT v_z AS A FUNCTION OF z WITH INTERACTION
plt.figure()
plt.title("With interaction")
plt.plot(r1RK_interacting[2,:], v1RK_interacting[2,:], label = "$P1$")
plt.plot(r2RK_interacting[2,:], v2RK_interacting[2,:], label = "$P2$")
plt.xlabel("$z [\mu m]$"); plt.ylabel("$v_z [\mu m/ \mu s]$")
plt.legend()
plt.show()

# PLOT 3D TRAJECTORY WITHOUT INTERACTION
fig = plt.figure()
ax1 = plt.axes(projection='3d')
ax1.set_title("Without interaction")
ax1.plot(r1RK[0,:], r1RK[1,:], r1RK[2,:], label='$P1$')
ax1.plot(r2RK[0,:], r2RK[1,:], r2RK[2,:], label='$P2$')
plt.legend()
plt.show()

# PLOT 3D TRAJECTORY WITH INTERACTION
fig = plt.figure()
plt.title("With interaction")
ax2 = plt.axes(projection='3d')
ax2.set_title("With interaction")
ax2.plot(r1RK_interacting[0,:], r1RK_interacting[1,:], r1RK_interacting[2,:], label='$P1$')
ax2.plot(r2RK_interacting[0,:], r2RK_interacting[1,:], r2RK_interacting[2,:], label='$P2$')
plt.legend()
plt.show()
"""

# PLOT RELATIVE ERROR BETWEEN ANALYTICAL AND NUMERICAL SOLUTION 
plt.figure()
plt.subplot(211)
plt.title('Error Euler')
plt.plot(t0, errEu0, label='n1 = 4000')
plt.plot(t1, errEu1, label='n2 = 8000')
plt.plot(t2, errEu2, label='n3 = 16000')
plt.plot(t3, errEu3, label='n4 = 32000')
plt.xlabel('$t [\mu s]$'); plt.ylabel('Relative error')
plt.legend()
plt.subplot(212)
plt.title('Error RK4')
plt.plot(t0, errRK0, label='n1 = 4000')
plt.plot(t1, errRK1, label='n2 = 8000')
plt.plot(t2, errRK2, label='n3 = 16000')
plt.plot(t3, errRK3, label='n4 = 32000')
plt.xlabel('$t [\mu s]$'); plt.ylabel('Relative error')
plt.legend()


"""
# PLOT PARTICLES LEFT
plt.figure()
plt.title('Particles left')
plt.plot(omV, NumParticles0/NumParticles0[0], 'o', label='f = 0.1')
plt.plot(omV, NumParticles1/NumParticles1[0], 'o', label='f = 0.4')
plt.plot(omV, NumParticles2/NumParticles2[0], 'o', label='f = 0.7')
plt.xlabel('$\omega_V [Hz]$'); plt.ylabel('Fraction of particles left')
plt.legend()
"""
plt.figure()
plt.plot(t0, r0[0])
plt.plot(t3, r1RK3[0,:])
plt.figure()
plt.plot(t0, r0[1])
plt.plot(t3, r1RK3[1,:])
plt.figure()
plt.plot(t0, r0[2])
plt.plot(t3, r1RK3[2,:])
