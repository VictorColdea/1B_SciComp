import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from datetime import datetime

class SolarSystem():
    """This class creates the SolarSystem object."""

    def __init__(self):
        """With self, you can access private attributes of the object."""
        self.size = 1000 # size of figure
        self.planets = []
        self.planet_lines = []
        # This initializes the 3D figure
        self.fig = plt.figure()
        self.ax = self.fig.subplot_mosaic([['Left','TopRight'],['Left','BottomRight']], gridspec_kw={'width_ratios':[2, 1]})
        #self.fig, self.ax = plt.subplots()
        #plt.subplot_mosaic([['left', 'right']], layout='constrained')
        self.ax_simulation = self.ax['Left']
        self.ax_momentum = self.ax['TopRight']
        self.ax_energy = self.ax['BottomRight']
        # array of times/frames
        self.time = 0

        # self.fig.add_axes(self.ax_simulation)
        self.dT = 1 # time increment

    def add_planet(self, planet):
        """Every time a planet is created, it gets put into the array."""
        self.planets.append(planet)
        momentum_line, = self.ax_momentum.plot([], [], lw=2)
        kinetic_energy_line, = self.ax_energy.plot([], [], lw=2, color='r')
        potential_energy_line, = self.ax_energy.plot([], [], lw=2, color='b')
        self.planet_lines.append({"momentum_line": momentum_line, "KE_line":kinetic_energy_line, "GPE_line":potential_energy_line})

    def update_planets(self):
        """This method moves and draws all of the planets."""
        self.ax_simulation.clear()
        self.time += 1
        #total_ang_mom, total_energy = 0, 0
        for i, planet in enumerate(self.planets[:-1]): # for each planet, carry out the move and draw commands
            planet.move()
            planet.draw()
            print(planet.kinetic_energies)
            assert len(planet.kinetic_energies) >= 1
            # planet.quantities_calc()
            planet.draw_quantities(**(self.planet_lines[i]))
        self.planets[-1].draw()
    

    # def update_plots(self):
    #     self.ax_momentum.plot(self.time, )
        
    ####    # for planet
        # total_ang_mom = self.planets[0].angular_momentum + self.planets[1].angular_momentum
        # total_energy = self.planets
        # print(f"Total ang. mom. = {total_ang_mom}")


    def fix_axes(self, frame_num):
        """The axes would change with each iteration otherwise."""
        # for simulation
        self.ax_simulation.set_xlim((-self.size/2, self.size/2)) # set the axes to be a certain size
        self.ax_simulation.set_ylim((-self.size/2, self.size/2))
        self.ax_simulation.set_xlabel("x")
        self.ax_simulation.set_ylabel("y")
        self.ax_simulation.set_title("Simulation of planets orbiting the Sun")
        self.ax_simulation.legend()
        for planet in range(len(self.planets)-1): #Omitting the sun (the last one)
            self.ax_simulation.text(self.planets[planet].position[0]-50,self.planets[planet].position[1]+20,f"Planet {planet+1}")
        self.ax_simulation.text(-25,40,"Sun")
        self.ax_simulation.text(400,-400,f"t={frame_num}")

        for ax in [self.ax_momentum, self.ax_energy]:
            # Things shared between the plots
            ax.set_xlim(0,frame_num)

        # Things not shared between the plots
        self.ax_momentum.set_title("Momentum plot")
        
        self.ax_energy.set_title("Energy plot")

        # for planet in self.planets[:-1]:
        #     planet.draw_quantities(momentum_line, kinetic_energy_line, potential_energy_line)
        #     #self.ax_momentum.plot(momentum_line)

                




    def gravity_planets(self):
        """This method calculates gravity interaction for every planet."""
        for i, first in enumerate(self.planets):
            for second in self.planets[i+1:]:
                first.gravity(second)

class Planet():
    """This class creates the Planet object."""

    def __init__(self, SolarSys, mass, num, position=(0, 0), velocity=(0, 0)): # starts off at origin at rest
        self.SolarSys = SolarSys # SolarSystem class
        self.mass = mass
        self.position = position
        self.velocity = velocity
        self.angular_momentums, self.kinetic_energies, self.potential_energies = [], [], []
        # The planet is automatically added to the SolarSys.
        self.SolarSys.add_planet(self)
        self.num = num
        #colors = ("green","blue")
        if self.num == 1:
            self.color = "green"
        else:
            self.color = "blue"
        

        #self.color = colors[1]

    def move(self):
        """The planet is moved based on the velocity."""
        self.position = (
            self.position[0] + self.velocity[0] * self.SolarSys.dT,
            self.position[1] + self.velocity[1] * self.SolarSys.dT
        )

    def draw(self):
        """The method to draw the planet."""
        # First to calculate the quantities
        self.quantities_calc()
        # Then, to draw the planet
        self.SolarSys.ax_simulation.plot(
            *self.position,
            marker="o",
            markersize=10,
            color=self.color,
            label=f"""L_{self.num} = {self.angular_momentum} \nKE_{self.num} = {self.kinetic_energy} \nGPE_{self.num} = {self.potential_energy}"""
        )
     
    def gravity(self, other):
        """The method to compute gravitational force for two planets. Numpy module is used to handle vectors."""
        distance = np.subtract(other.position, self.position) # define r as displacement from self to other
        distanceMag = np.linalg.norm(distance) # calculate magnitude of distance
        distanceUnit = np.divide(distance, distanceMag) # unit vector pointing in right direction
        forceMag = self.mass * other.mass / (distanceMag ** 2) # gravity force formula (magnitude)
        force = np.multiply(distanceUnit, forceMag) # force vector

        # Switch makes force on self opposite to other
        switch = 1
        for body in (self, other):
            acceleration = np.divide(force, body.mass) #this is accel. due to the grav. force
            velocity_change = np.multiply(acceleration, self.SolarSys.dT * switch)
            body.velocity = np.add(body.velocity, velocity_change)
            switch *= -1 # so it's different to the "other" planet

    def angular_momentum_calc(self):
        # Calculating angular momentum wrt Sun. (ignoring with other planets)
        # w = (r x v)
        # cross product comes from r1v2-r2v1
        ang_vel_1 = np.multiply(self.position[0],self.velocity[1])
        ang_vel_2 = np.multiply(self.position[1],self.velocity[0])
        ang_velocity = np.subtract(ang_vel_1,ang_vel_2)

        # L = m w
        self.angular_momentum = np.multiply(ang_velocity, self.mass)
        return self.angular_momentum

    def kinetic_energy_calc(self):
        # Calculating kinetic energy
        # KE = 1/2 * m * v**2
        velocityMag = np.linalg.norm(self.velocity)
        self.kinetic_energy = 1/2 * np.multiply(self.mass, velocityMag**2)
        return self.kinetic_energy

    def potential_energy_calc(self):
        # Calculating potential energy - this is odd because we have neglected G
        # GPE = - GMm/r traditionally.
        # However, we can just take GPE = -Mm/r since G = 1 in this simulation
        massMultiply = self.mass * 10000
        distance_from_sun = self.position # define r as displacement from self to other
        distanceMag = np.linalg.norm(distance_from_sun) # calculate magnitude of distance
        self.potential_energy = - massMultiply / distanceMag
        return self.potential_energy

    def quantities_calc(self):
        #calculate quantities of the planet for future plots
        self.angular_momentums.append(self.angular_momentum_calc())
        # print(self.kinetic_energies)
        self.kinetic_energies.append(self.kinetic_energy_calc())
        # print(self.kinetic_energies)
        self.potential_energies.append(self.potential_energy_calc())

    def draw_quantities(self, momentum_line, KE_line, GPE_line):
        time_data = range(self.SolarSys.time)
        assert len(self.angular_momentums) == len(time_data), f"len quantity: {len(self.angular_momentums)}, len time: {len(time_data)}"
        momentum_line.set_data(time_data, self.angular_momentums)
        KE_line.set_data(time_data, self.kinetic_energies)
        GPE_line.set_data(time_data, self.potential_energies)




class Sun(Planet):
    """This class is inherited from Planet. Everything is the same as in the 
    planet, except that the position of the sun is fixed. Also, the color is yellow."""
    def __init__(self, SolarSys, mass=10000, position=(0, 0), velocity=(0, 0)):
        super(Sun, self).__init__(SolarSys, mass, position, velocity)
        self.color = "orange"
        self.num = "Sun"

    def move(self):
        self.position = self.position # overrides the planet.move command so that the sun doesn't move

    def draw(self):
        """The method to draw the Sun."""
        # Then, to draw the planet
        self.SolarSys.ax_simulation.plot(
            *self.position,
            marker="o",
            markersize=20,
            color=self.color,
        )

# Instantiating of the solar system.
SolarSys = SolarSystem()

# Instantiating of planets.

planet1 = Planet(SolarSys, mass=10, num=1, position=(100, 0), velocity=(0, 10))
planet2 = Planet(SolarSys, mass=10, num=2, position=(-200, 0), velocity=(0, -8))
# planet3 = Planet(SolarSys, mass=10, position=(-100, -50, 150), velocity=(-4, 0, 0))
# planet4 = Planet(SolarSys, mass=10, position=(-200, -50, 150), velocity=(-2, 2, 0))

# Instantiating of the sun.
sun = Sun(SolarSys)

def animate(i):
    """This controls the animation."""
    print("The frame is:",i)
    SolarSys.gravity_planets()
    SolarSys.update_planets()
    SolarSys.fix_axes(i)


# This calls the animate function and creates an animation.
anim = animation.FuncAnimation(SolarSys.fig, animate, frames=100, interval=100)

plt.show()

date = datetime.now()

# This prepares the writer for the animation.
#writervideo = animation.FFMpegWriter(fps=60)

# This saves the animation.
# time = datetime.now().strftime("%X")[:5]
# #file_loc = "C:\Users\victo\OneDrive - University of Cambridge\Year2\Computing\SciComp\Session5_Planets\Planets_attempts\planets_animation"
# anim.save(f"animation_.gif", writer="pillow", dpi=200)
