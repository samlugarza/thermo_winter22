'''
This code is to solve question 2 (a and b) for HW 1 of Thermodynamics Winter 2022
'''

# import the code for the simulation
from simulator import Simulation
import matplotlib.pyplot as plt
import numpy as np
from time import time

# create a new class for the simulation with some randomish variable choices
#Set visualse = True to see the simulation for part a; false to seee graphs faster in part b 

masses = np.ones((1,100))

sim = Simulation(N=100, E=1, size=750, masses = masses[0], rad=3, delay=5, visualise=False) 

'''
Units using for the simulation above: for reference and erg = cm^2 * g * s^-2

N is unitless; E = [erg]; size = [cm]; mass = [g]; rad = [cm]; and the delay is milliseconds [ms]
'''

analytic_derivation = True  #means that I still have to do the math for this later after I finish the code 

analytic_implementation = True # But they looky wonky and need some definite work 

start = time()

'''
This code is to solve question 2 part b 
'''

# run the simulation for a long time step 
sim.run_simulation(steps=10000) 

velocity_x = [sim.particles[i].vx for i in range(sim.N)]
velocity_y = [sim.particles[i].vy for i in range(sim.N)]

#Then want to run the simulation for a lot of particles for a short time while saving the velocities for each and then combine so that have nice histogram 
for i in range(100): 
    sim.run_simulation(steps=100)
    x_velocity = [sim.particles[i].vx for i in range(sim.N)]
    y_velocity = [sim.particles[i].vy for i in range(sim.N)]

    velocity_x.extend(x_velocity)
    velocity_y.extend(y_velocity)

print("Runtime: {:1.2f}".format(time() - start))

velocity_x = np.array(velocity_x)
velocity_y = np.array(velocity_y)

speed = np.sqrt(velocity_x**2 + velocity_y**2)

print('speed', speed)
print(speed.shape)

np.save('speeds.npy', speed) #this saves the speeds and then I make the plots in a jupyter notebook since it is not plotting here 

#Now I need to find the energy in ergs so that I can plot that as a histogram 

masses = np.ones((1,10100))
E = (masses[0] * speed**2)/2 

print('Energy', E)

np.save('energy.npy', E) #this saves the energies and then I am making the plots in a jupyter notebook 


