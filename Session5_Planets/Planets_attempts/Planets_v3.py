import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

class SolarSystem():
    """This class creates the SolarSystem object."""

    def __init__(self):
        """With self, you can access private attributes of the object."""
        self.size = 1000 # size of figure
        self.planets = []
        # This initializes the 3D figure
        self.fig, self.ax = plt.subplots()
        #plt.subplot_mosaic([['left', 'right']], layout='constrained')
        self.fig.add_axes(self.ax)
        self.dT = 1 # time increment

    def add_planet(self, planet):
        """Every time a planet is created, it gets put into the array."""
        self.planets.append(planet)

    def update_planets(self):
        """This method moves and draws all of the planets."""
        self.ax.clear()
        total_ang_mom, total_energy = 0, 0
        for planet in self.planets: # for each planet, carry out the move and draw commands
            planet.move()
            planet.draw()


        
        #for planet
        # total_ang_mom = self.planets[0].angular_momentum + self.planets[1].angular_momentum
        # total_energy = self.planets
        # print(f"Total ang. mom. = {total_ang_mom}")


    def fix_axes(self):
        """The axes would change with each iteration otherwise."""
        #self.fig.add_subplot(1,2,1)
        self.ax.set_xlim((-self.size/2, self.size/2)) # set the axes to be a certain size
        self.ax.set_ylim((-self.size/2, self.size/2))
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.set_title("Simulation of planet orbiting the Sun")
        self.ax.legend()
        for planet in range(len(self.planets)-1): #Omitting the sun (the last one)
            self.ax.text(self.planets[planet].position[0]-50,self.planets[planet].position[1]+20,f"Planet {planet+1}")
        self.ax.text(-25,20,"Sun")


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
        self.angular_momentum_calc()
        self.kinetic_energy_calc()
        self.potential_energy_calc()
        # Then, to draw the planet
        self.SolarSys.ax.plot(
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

    def kinetic_energy_calc(self):
        # Calculating kinetic energy
        # KE = 1/2 * m * v**2
        velocityMag = np.linalg.norm(self.velocity)
        self.kinetic_energy = 1/2 * np.multiply(self.mass, velocityMag**2)

    def potential_energy_calc(self):
        # Calculating potential energy - this is odd because we have neglected G
        # GPE = - GMm/r traditionally.
        # However, we can just take GPE = -Mm/r since G = 1 in this simulation
        massMultiply = self.mass * 10000
        distance_from_sun = self.position # define r as displacement from self to other
        distanceMag = np.linalg.norm(distance_from_sun) # calculate magnitude of distance
        self.potential_energy = - massMultiply / distanceMag






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
        self.SolarSys.ax.plot(
            *self.position,
            marker="o",
            markersize=10,
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
    SolarSys.fix_axes()

# This calls the animate function and creates an animation.
anim = animation.FuncAnimation(SolarSys.fig, animate, frames=100, interval=100)
plt.show()

# This prepares the writer for the animation.
#writervideo = animation.FFMpegWriter(fps=60)

# This saves the animation.
#anim.save("planets_animation.gif", writer="pillow", dpi=200)
