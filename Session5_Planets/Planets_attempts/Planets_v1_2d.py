import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
#from mpl_toolkits.mplot3d import Axes3D

class SolarSystem():
    """This class creates the SolarSystem object."""

    def __init__(self):
        """With self, you can access private attributes of the object."""
        self.size = 1000 # size of figure
        self.planets = []
        # This initializes the 3D figure
        self.fig, self.ax = plt.subplots()
        # self.fig = plt.figure()
        # self.ax = Axes3D(self.fig, auto_add_to_figure=False)
        self.fig.add_axes(self.ax)
        self.dT = 1 # time increment

    def add_planet(self, planet):
        """Every time a planet is created, it gets put into the array."""
        self.planets.append(planet)

    def update_planets(self):
        """This method moves and draws all of the planets."""
        self.ax.clear()
        for planet in self.planets: # for each planet, carry out the move and draw commands
            planet.move()
            planet.draw()

    def fix_axes(self):
        """The axes would change with each iteration otherwise."""
        self.ax.set_xlim((-self.size/2, self.size/2)) # set the axes to be a certain size
        self.ax.set_ylim((-self.size/2, self.size/2))
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.set_title("Simulation of planet orbiting the Sun")


    def gravity_planets(self):
        """This method calculates gravity interaction for every planet."""
        for i, first in enumerate(self.planets):
            for second in self.planets[i+1:]:
                first.gravity(second)

class Planet():
    """This class creates the Planet object."""

    def __init__(self, SolarSys, mass, position=(0, 0), velocity=(0, 0)): # starts off at origin at rest
        self.SolarSys = SolarSys # SolarSystem class
        self.mass = mass
        self.position = position
        self.velocity = velocity
        # The planet is automatically added to the SolarSys.
        self.SolarSys.add_planet(self)
        self.color = "black"

    def move(self):
        """The planet is moved based on the velocity."""
        self.position = (
            self.position[0] + self.velocity[0] * self.SolarSys.dT,
            self.position[1] + self.velocity[1] * self.SolarSys.dT#,
            #self.position[2] + self.velocity[2] * self.SolarSys.dT
        )

    def draw(self):
        """The method to draw the planet."""
        self.SolarSys.ax.plot(
            *self.position,
            marker="o",
            markersize=10,
            color=self.color
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
            acceleration = np.multiply(force, self.SolarSys.dT * switch)
            body.velocity = np.add(body.velocity, acceleration)
            switch *= -1 # so it's different to the "other" planet


class Sun(Planet):
    """This class is inherited from Planet. Everything is the same as in the 
    planet, except that the position of the sun is fixed. Also, the color is yellow."""
    def __init__(self, SolarSys, mass=1000, position=(0, 0), velocity=(0, 0)):
        super(Sun, self).__init__(SolarSys, mass, position, velocity)
        self.color = "yellow"

    def move(self):
        self.position = self.position # overrides the planet.move command so that the sun doesn't move

# Instantiating of the solar system.
SolarSys = SolarSystem()

# Instantiating of planets.

planet1 = Planet(SolarSys, mass=10, position=(150, 50), velocity=(-5, 5))
# planet2 = Planet(SolarSys, mass=10, position=(100, -50, 150), velocity=(5, 0, 0))
# planet3 = Planet(SolarSys, mass=10, position=(-100, -50, 150), velocity=(-4, 0, 0))
# planet4 = Planet(SolarSys, mass=10, position=(-200, -50, 150), velocity=(-2, 2, 0))

# Instantiating of the sun.
sun = Sun(SolarSys)

def animate(i):
    """This controls the animation."""
    print("The frame is:",i)
    SolarSys.gravity_planets()
    SolarSys.update_planets()
    SolarSys.fix_axes()

# This calls the animate function and creates an animation.
anim = animation.FuncAnimation(SolarSys.fig, animate, frames=100, interval=100)
plt.show()

# This prepares the writer for the animation.
#writervideo = animation.FFMpegWriter(fps=60)

# This saves the animation.
#anim.save("planets_animation.mp4", writer=writervideo, dpi=200)
