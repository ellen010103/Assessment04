import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)

#define lennard jones parameters for argon
kB = 1.38e-23 #J/K
eps = 125.7 * kB #J
sigma = 3.345e-10 #m

#define the timestep
dt = 1e-25

#define the size of the simulation box
box_size = 5*sigma   #m

#define the number of particles 
N = 2

#define the mass of the particles
mass = 6.6335209e-26 #kg, mass of argon atom

#initial conditions
pos = np.array([[0.5*box_size, 0.5*box_size],
              [0.5*box_size+sigma, 0.5*box_size+sigma]])   #initial positions, particle i located in centre of box and j located a disantce of sigma away
vel = np.random.rand(N,2) #initial velocities

#define a function to calculate the lennard jones potential
def LJ_potential(pos):                     
    return 4*eps*((sigma/pos)**12 - (sigma/pos)**6)   

#define a function to calculate the force between two particles using the differential of the lennard jones potential
def LJ_force(pos):                      
    return 24*(eps/pos)*(-2*((sigma/pos)**12) + (sigma/pos)**6)   

#define a function to calculate the net force on each particle
def calculate_forces(pos):
    F_net = np.zeros((N, 2))           #empty array for adding calculate net force to
    for i in range(N):
        for j in range(i+1, N): 
            if i != j:
                diff_pos = pos[j] - pos[i]
                pos_ij = np.sqrt(np.sum(diff_pos**2)) #find distance between particles
                F_ij = LJ_force(pos_ij) * diff_pos/pos_ij
                F_net[i] += F_ij
                F_net[j] -= F_ij
    return F_net 

 
#define a function to perform the Verlet integration scheme
def verlet_integrate(pos, vel, dt, mass):
    acc = calculate_forces(pos) / mass
    pos_new = pos + dt*vel + 0.5 * acc*dt*dt
    if pos_new.any() <= 0.5*box_size:
        acc_new = (1.0/mass) * acc
        vel_new = vel + (dt/2) * (acc+acc_new)
    else:
        pos_new[0] = pos[0]*-1
    return pos_new, vel_new

#perform the Verlet integration for a specified number of timesteps
num_steps = 10000
position = np.zeros(num_steps)
time = np.zeros(num_steps)

for i in range(num_steps):                        #verlet_integrate function not used here to allow for if statement to work
    acc = calculate_forces(pos) / mass
    pos_new = pos + dt*vel + 0.5 * acc*dt*dt
    if pos_new.any() <= 0.5*box_size:             #if particle is within box continue
        acc_new = (1.0/mass) * acc
        vel_new = vel + (dt/2) * (acc+acc_new)

    position[i]=((np.sqrt(np.sum((pos_new[0]-pos_new[1])**2))))
    time[i] = i * dt
        
position = position/sigma
print(position)

import matplotlib.pyplot as plt
#plot against time to see how the distnace beeen the particles changes over time
fig, ax = plt.subplots()
ax.plot(time, position)

#plot distance between particles against lennard jones potential
fig, ax = plt.subplots()
ax.plot(position, LJ_potential(position))