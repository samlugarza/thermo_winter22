from turtle import position
import numpy as np
import tkinter as tk           # simple gui package for python
import scipy 
import scipy.stats

class particle():
    def __init__(self, size, pid, mass, init_posx=None, init_posy=None, init_v=5, rad=3):
        """Initialise the particles

        Parameters
        ----------
        size : int
            Size of the box
        pid : int
            Unique particle ID
        init_v : int, optional
            Initial velocity for particles, by default 5
        rad : int, optional
            Radius of each particle, by default 3

        mass: 
            mass of the particle
        """

        kB = 1.38e-16     # Boltzmann const. [cm^2 g s^-2 K^-1]

        # choose random x and y positions within the grid (padded by radius of particles)
        self.x = init_posx #np.random.uniform(0 + rad, size - rad)
        self.y = init_posy #np.random.uniform(0 + rad, size - rad)

        # set random velocities for each particle (randomly distributed between x and y speed)
        self.vx = np.random.uniform(0, init_v) * np.random.choice([-1, 1])
        self.vy = np.sqrt(init_v**2 - self.vx**2) * np.random.choice([-1, 1])

        # set the radius of the particle
        self.rad = rad

        # assign a particle id to each particle
        self.pid = pid

        # Set the mass of the particle
        self.mass = mass 

    def update_x(self, val):
        self.x = val

    def update_y(self, val):
        self.y = val

    def update_vx(self, val):
        self.vx = val

    def update_vy(self, val):
        self.vy = val


class Simulation():  # this is where we will make them interact
    def __init__(self, N, E, size, masses, rad, init_posx=None, init_posy=None, delay=20, visualise=True):
        """Simulation class initialisation. This class handles the entire particle
        in a box thing.

        Parameters
        ----------
        N : `int`
            Total number of particles
        E : `int`
            Kinetic energy to start with
        size : `int`
            Size of the box
        rad : `int`
            Radius of the particles

        mass: 
            mass of the particle 
            
        delay : `int`
            Delay in milliseconds between showing/running timesteps
        """
        
        # Physical contants
        kB = 1.38e-16   # Boltzmann const. [cm^2 g s^-2 K^-1]
        mH = 1.67e-24   # Atomic Hydrogen mass [g]
        
        self.N = N
        self.E = E
        self.size = size
        self.rad = rad
        self.delay = delay 

        # create an empty position array and try to find a position for each particle
        positions = []
        for _ in range(N):
            # keep looping until the position is no longer invalid
            position_invalid = True
            while position_invalid:
                # pick a random (x, y) position
                possible_pos = np.random.uniform(0 + rad, size - rad, size=2)
                position_invalid = False

                # loop over all other chosen positions
                for other in positions:
                    # mark the position as bad if it overlaps with another particle
                    position_invalid = np.sqrt(sum((other - possible_pos)**2)) <= 2 * rad
                    if position_invalid:
                        break

            # add to the position array
            positions.append(possible_pos)

        # initialise N particle classes
        self.particles = [particle(size=size, pid=i, mass=masses[i], init_v=E, rad=rad, init_posx=positions[i][0], init_posy=positions[i][1]) for i in range(N)]

        self.visualise = visualise
        
        if visualise:
            self.delay = delay
            self.canvas = None
            self.root = None
            self.particle_handles = {}

            self._init_visualization()
            self.root.update()

    def _init_visualization(self):
        # start the visualisation box
        self.root = tk.Tk()
        self.root.title("Particles in a Box!")

        # create a canvas with the right size
        self.canvas = tk.Canvas(self.root, width=self.size, height=self.size)
        self.canvas.pack()

        # add a close button
        self.button = tk.Button(self.root, text='Close', command=self._quit_visualisation)
        self.button.place(x=self.size, y=10, anchor="e")

        self.timestep_message = self.canvas.create_text(self.size // 2, 10, text="Timestep = 0")

        # add all of the particles
        for p in self.particles:
            self.particle_handles[p.pid] = self._draw_particle(p)

        # update this all on the canvas
        self.root.update()

    def _quit_visualisation(self):
        self.root.destroy()

    def _draw_particle(self, particle):
        """Draw a circle on the canvas corresponding to particle

        Returns the handle of the tkinter circle element"""
        x0 = particle.x - particle.rad
        y0 = particle.y - particle.rad
        x1 = particle.x + particle.rad
        y1 = particle.y + particle.rad
        return self.canvas.create_oval(x0, y0, x1, y1, fill='black', outline='black')

    def _move_particle(self, particle):
        xx = particle.x + particle.vx
        yy = particle.y + particle.vy
        particle.update_x(xx)
        particle.update_y(yy)

        if self.visualise:
            self.canvas.move(self.particle_handles[particle.pid], particle.vx, particle.vy)

    def resolve_particle_collisions(self):
        #raise NotImplementedError
        for i in range(self.N):
            for j in range(i + 1, self.N):
                distance = self._calculate_distance(self.particles[i], self.particles[j]) 
                if distance <= self.particles[i].rad + self.particles[j].rad:
                    # do stuff and there was much rejoicing.... yay :D
                    #These are the 2D collision equations with two moving objects from: 
                    #   https://en.wikipedia.org/wiki/Elastic_collision#Two-dimensional_collision_with_two_moving_objects

                    x1 = np.array([self.particles[i].x, self.particles[i].y])
                    x2 = np.array([self.particles[j].x, self.particles[j].y])

                    v1 = np.array([self.particles[i].vx, self.particles[i].vy])
                    v2 = np.array([self.particles[j].vx, self.particles[j].vy])

                    v1_f = v1 - (2 * self.particles[j].mass / (self.particles[i].mass + self.particles[j].mass)) *\
                        np.dot(v1 - v2, x1 - x2) / np.linalg.norm(x1 - x2)**2 * (x1 - x2)
                    
                    v2_f = v2 - (2 * self.particles[i].mass / (self.particles[i].mass + self.particles[j].mass)) *\
                        np.dot(v2 - v1, x2 - x1) / np.linalg.norm(x2 - x1)**2 * (x2 - x1)
                    
                    self.particles[i].update_vx(v1_f[0])
                    self.particles[i].update_vy(v1_f[1])

                    self.particles[j].update_vx(v2_f[0])
                    self.particles[j].update_vy(v2_f[1]) 


    def _calculate_distance(self, particle1, particle2):
        dist = np.sqrt((particle1.x - particle2.x)**2 + (particle1.y - particle2.y)**2)
        return dist


    def resolve_wall_collisions(self):
        # check whether each particle hits the wall
        # for each collider reflect its velocity (account for ones that hit both walls)
        #raise NotImplementedError

        for particle in self.particles: 
            if particle.x >= self.size - particle.rad or particle.x <= 0 + particle.rad: 
                particle.update_vx(-particle.vx)

            if particle.y >=self.size - particle.rad or particle.y <= 0 + particle.rad:
                particle.update_vy(-particle.vy)
    
    '''
    #This is the code that I was working on for Question 3 

    steps = []
    
    def _cdf(self):
        #calculate the velocties (x and y) for each of the particles 
        velocity_x = [self.particles[i].vx for i in range(self.N)]
        velocity_y = [self.particles[i].vy for i in range(self.N)]
        #put those values into arrays so that I can work with them 
        velocity_x = np.array(velocity_x)
        velocity_y = np.array(velocity_y)
        #calculate the speed 
        speed = np.sqrt(velocity_x**2 + velocity_y**2)
        #need to find kT in terms of the vrms because don't have a set temperature 
        vrms = np.sqrt(np.mean(speed)**2)
        k_BT = (self.particles.mass*vrms**2)/3 
        
        cdf_dNdv = 1-np.exp(-(self.particles.mass/(2*k_BT))*speed**2)
        #cdf_dNdE = 1-np.exp(-self.particles.E/k_BT)
        
        return cdf_dNdv, #cdf_dNdE

    def _ks_test(self): 
        #calculate the velocties (x and y) for each of the particles 
        velocity_x = [self.particles[i].vx for i in range(self.N)]
        velocity_y = [self.particles[i].vy for i in range(self.N)]
        #put those values into arrays so that I can work with them 
        velocity_x = np.array(velocity_x)
        velocity_y = np.array(velocity_y)
        #calculate the speed 
        speed = np.sqrt(velocity_x**2 + velocity_y**2)
        #need to find kT in terms of the vrms because don't have a set temperature 
        vrms = np.sqrt(np.mean(speed)**2)
        k_BT = (self.particles.mass*vrms**2)/3 

        #need to run a KS test of the speed and the energy 
        kstest_speed, pval = scipy.stats.kstest(speed, self.cdf_dNdv)
        #kstest_energy = scipy.stats.kstest()

        return kstest_speed, pval

    def steady_state(self): 
        pthreshold = 0.05 
        
        for i in range(steps)
        
        if pval > pthreshold: 
            #need ti find a way to add a value to the existing list of steps and keep going 

        if pval <= pthreshold: 
            #this is where the appending to the time step would end and the simualtion would stop running 

    '''


    def run_simulation(self, steps=1000):
        for i in range(steps):
            # 1. update all particle positions based on current speeds
            for particle in self.particles:
                self._move_particle(particle)

            # 2. resolve whether any hit the wall and reflect them
            self.resolve_wall_collisions()

            # 3. resolve any particle collisions and transfer momentum
            self.resolve_particle_collisions()

            # update visualization with a delay
            if self.visualise: #this is an idea from Tom so that I can plot the graphs faster in problem_1.py
                self.root.after(self.delay, self.root.update())

                # change the timestep message as well
                self.canvas.itemconfig(self.timestep_message, text="Timestep = {}".format(i))

        if self.visualise: #this is an idea from Tom so that I can plot the graphs faster in problem_1.py
            self.root.mainloop()

    def get_velocities(self):
        raise NotImplementedError
