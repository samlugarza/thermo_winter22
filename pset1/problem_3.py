'''
This code is to solve question 3 (a and b and c) for HW 1 of Thermodynamics Winter 2022

Relaxation Timescale: Initially the distribution of particle energies changes rapidly, but eventually it approaches a steady-state.

'''

# import the code for the simulation
from simulator import Simulation
import matplotlib.pyplot as plt
import numpy as np
from time import time
from scipy import optimize
from scipy.stats import cauchy, norm, sem
import scipy.stats

# create a new class for the simulation with some randomish variable choices
#Set visualse = True to see the simulation for part a; false to seee graphs faster in part b 

masses = np.ones((1,100))

sim = Simulation(N=100, E=1, size=750, masses = masses[0], rad=3, delay=5, visualise=False) 

'''
Units using for the simulation above: for reference and erg = cm^2 * g * s^-2

N is unitless; E = [erg]; size = [cm]; mass = [g]; rad = [cm]; and the delay is milliseconds [ms]
'''

start = time()

# run the simulation for a long time step 
sim.run_simulation(steps=100) 

#Just trying some stuff out and want to see if I can make the graphs for part b 
#assuming the cross section is the radius? Since it is units of length???? 

C = 1.2
velocity_x = np.array([sim.particles[i].vx for i in range(sim.N)])
velocity_y = np.array([sim.particles[i].vy for i in range(sim.N)])

speed = np.sqrt(velocity_x**2 + velocity_y**2)
print('speed', speed.shape)

cdfspeed = speed.cdf(0)

print('speedcdf', cdfspeed)

'''
trelax = C/(sim.N*sim.rad*speed)
print('trelax', trelax.shape)

N = np.array(range(1,101))
print('N', N)

plt.plot(N,trelax)
plt.show()'''