'''
This code is to solve question 1 (a and b) for HW 1 of Thermodynamics Winter 2022
'''
# import the code for the simulation
from simulator import Simulation
import matplotlib.pyplot as plt
import numpy as np
from time import time

# create a new class for the simulation with some randomish variable choices
#Set visualse = True to see the simulation for part a; false to seee graphs faster in part b 

masses = np.ones((1,100))*4

sim = Simulation(N=100, E=3, size=750, masses = masses[0], rad=15, delay=5, visualise=True) 

'''
Units using for the simulation above: for reference and erg = cm^2 * g * s^-2

N is unitless; E = [erg]; size = [cm]; mass = [g]; rad = [cm]; and the delay is milliseconds [ms]
'''

#Use the bottom code to aswer part b 
#1) Intital Position Scatter Plot 
positions_x = [sim.particles[i].x for i in range(sim.N)]
positions_y = [sim.particles[i].y for i in range(sim.N)]

fig, axes = plt.subplots(2, 2, figsize=(10, 10))

axes[0, 0].scatter(positions_x, positions_y)
axes[0, 0].set_xlabel('x [cm]')
axes[0, 0].set_ylabel('y [cm]')
axes[0, 0].set_title('Initial Position')

 #2) Initial Speed Histogram 
 #Need to turn them into a speed? and not have the independent x and y but the total one and then make a speed. 

velocity_x = np.array([sim.particles[i].vx for i in range(sim.N)])
velocity_y = np.array([sim.particles[i].vy for i in range(sim.N)])

speed = np.sqrt(velocity_x**2 + velocity_y**2)

axes[0, 1].hist(speed)
axes[0, 1].set_xlabel('Speed, v [cm/s]')
axes[0, 1].set_title('Initial Speed')


# run the simulation
sim.run_simulation(steps=100) 

#3) Final Position Scatter Plot 
positions_x = [sim.particles[i].x for i in range(sim.N)]
positions_y = [sim.particles[i].y for i in range(sim.N)]

axes[1, 0].scatter(positions_x, positions_y)
axes[1, 0].set_xlabel('x [cm]')
axes[1, 0].set_ylabel('y [cm]')
axes[1, 0].set_title('Final Position')


#4) Final Speed Histogram 
velocity_x = np.array([sim.particles[i].vx for i in range(sim.N)])
velocity_y = np.array([sim.particles[i].vy for i in range(sim.N)])

speed = np.sqrt(velocity_x**2 + velocity_y**2)

axes[1, 1].hist(speed)
axes[1, 1].set_xlabel('Speed, v [cm/s]')
axes[1, 1].set_title('Final Speed')
plt.show()






