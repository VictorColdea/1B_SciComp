import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import itertools

class Canvas():
    """This class creates the canvas object."""

    def __init__(self):
        """With self, you can access private attributes of the object."""
        self.size = 20
        self.balls = []
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot()
        self.dT = 1 #time of each frame - making this smaller is EXACTLY the same as making the speed smaller

    def add_ball(self, ball):
        """Every time a ball is created, it gets put into the array."""
        self.balls.append(ball)

    def update_balls(self):
        """This method moves and draws all balls."""
        self.ax.clear()
        for i, ball in enumerate(self.balls):
            ball.move()
            ball.draw()
            
            

    def fix_axes(self):
        """The axes would change with each iteration otherwise."""
        self.ax.set_xlim((-self.size/2, self.size/2))
        self.ax.set_ylim((-1, 1))
        self.ax.set_title("Simulation of particles confined in 1D")
        # for i, ball in enumerate(self.balls):
        #     self.ax.text(ball.position, 10, f"Ball {i}")   ####why doeS THIS NOT WORK

    def check_collision(self):
        """This method checks if balls are colliding."""
        combinations = list(itertools.combinations(range(len(self.balls)), 2)) #creates a list of all the combinations between balls
        for pair in combinations:
            self.balls[pair[0]].collide(self.balls[pair[1]])

        #also need to consider collisions between balls and "walls"

class Ball():
    """This class creates the ball object."""
    #should be ball icl

    def __init__(self, canvas, mass, position=0, velocity=0):
        self.canvas = canvas
        self.mass = mass
        self.size = 10*np.sqrt(self.mass)
        self.position = position
        self.velocity = velocity
        # The ball is automatically added to the canvas.
        self.canvas.add_ball(self)
        #self.color = "black"

        #Defining radius
        r0, a0 = 2.5, 90
        const_of_prop = r0/a0 #constant of proportionality between size and radius
        self.radius = self.size * const_of_prop

    def move(self):
        """The ball is moved based on the velocity."""
        self.position = self.position + self.velocity * self.canvas.dT
        self.kinetic_energy = 1/2 * self.mass * (self.velocity)**2

    def draw(self):
        """The method to draw the ball. Note: if you don’t
        specify the color of the ball, the color of each new
        element that is plotted will be randomly assigned.
        """
        self.canvas.ax.plot(self.position, 0,"o", markersize=self.size)#, c=self.color) # moves with xpos = self.position and ypos = 0 (i.e. moves in 1D)
        

    def collide(self, other):
        # Watch out with the threshold, this should
        # be greater than velocity. If, in a single
        # iteration, balls move relatively more than
        # what the threshold is, the collision might
        # not be triggered.

        #### we can define the threshold depending on the speed and radius of the particles
        centres_dist = self.radius + other.radius
        #speed_of_approach = abs(self.velocity-other.velocity)
        #so in one frame, centres_dist changes by -speed_of_appr * dT (if they are moving together)


          
        threshold = centres_dist# - speed_of_approach*self.canvas.dT
        #^ the effect of speed_of_approach on the threshold I am unsure of.

        if abs(self.position - other.position) < threshold:
            # The velocity after collision is not modeled
            # correctly. Take a look at what happens when
            # particles going in the same direction collide.
            # How to fix that? 
            # ^ good question matey
            u1, u2 = self.velocity, other.velocity
            m1, m2 = self.mass, other.mass
            v1 = ((1-m2/m1)/(1+m2/m1))*u1 + (2/(1+m1/m2))*u2
            v2 = (2/(1+m2/m1))*u1 - ((1-m2/m1)/(1+m2/m1))*u2

            
            self.velocity = v1
            other.velocity = v2

        # Wall collision
        for ball in [self, other]:
            if abs(ball.position) > self.canvas.size/2 - ball.radius:
                ball.velocity *= -1 #bounce off the wall if radius touches it

canvas = Canvas()
ball1 = Ball(canvas, mass=1, position=-2, velocity=0.1)
ball2 = Ball(canvas, mass=4, position=2, velocity=-0.2)
ball3 = Ball(canvas, mass=5, position=4, velocity=-0.03)
ball4 = Ball(canvas, mass=2, position=-5, velocity=0.05)
#momentums above around 1 cause issues

def animate(i):
    """This controls the animation."""
    #print("The frame is:", i)
    canvas.update_balls()
    canvas.check_collision()
    canvas.fix_axes()

# This calls the animate function and creates animation.
anim = animation.FuncAnimation(canvas.fig, animate, frames=500, interval=10)
plt.show()

# # This prepares the writer for the animation.
# writervideo = animation.FFMpegWriter(fps=60)

# # This saves the animation.
# anim.save("balls_animation.mp4", writer=writervideo, dpi=200)
