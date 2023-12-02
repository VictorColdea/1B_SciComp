import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import itertools

"""Written by Victor Coldea (vic21), 29/11/2023"""

class Canvas():
    """This class creates the canvas object."""

    def __init__(self):
        """With self, you can access private attributes of the object."""
        self.size = 20
        self.balls = []
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot()
        self.dT = 3.5  #time of each frame - making this smaller is EXACTLY the same as making the speed smaller
        self.pressure = 0


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
        self.ax.set_ylim((-self.size/2, self.size/2))
        self.ax.set_title("Simulation of particles confined in a 2D box")

    def check_collision(self):
        """This method checks if balls are colliding."""
        combinations = list(itertools.combinations(range(len(self.balls)), 2)) #creates a list of all the combinations between balls
        for pair in combinations:
            self.balls[pair[0]].collide(self.balls[pair[1]])

    def total_kinetic_energy_calc(self):
        total_kinetic_energy = 0
        for ball in self.balls:
            KE = ball.kinetic_energy_calc()
            total_kinetic_energy += KE
        self.total_KE = total_kinetic_energy

    # def total_momentum_calc(self):
    #     total_momentum = np.array([0,0])
    #     for ball in self.balls:
    #         P = ball.momentum_calc() #define this
    #         total_momentum = total_momentum + np.array(P)
    #     self.total_momentum = total_momentum

            # ^ this is meaningless due to the walls messing up conservation of momentum - at least to the computer  

class Ball():
    """This class creates the ball object."""
    #should be ball icl

    def __init__(self, canvas, mass, position=[0, 0], velocity=[0, 0]):
        self.canvas = canvas
        self.mass = mass
        self.size = 10*np.sqrt(self.mass)
        self.position = np.array(position)
        self.velocity = np.array(velocity)
        # The ball is automatically added to the canvas.
        self.canvas.add_ball(self)
        #self.color = "black"

        #Defining radius
        r0, a0 = 2.5, 90
        indice = 1
        const_of_prop = r0/(a0**indice) #constant of proportionality between  matplotlib size and radius
        self.radius = (self.size)**indice * const_of_prop
        

    def move(self):
        """The ball is moved based on the velocity."""
        self.position = self.position + self.velocity * self.canvas.dT

    def kinetic_energy_calc(self):
        """Calc. of KE"""
        #self.kinetic_energy
        KE = 1/2 * self.mass * (np.linalg.norm(self.velocity))**2
        return KE
    
    def momentum_calc(self):
        """Calc. of linear momentum"""
        P = self.mass * self.velocity
        return P

    def draw(self):
        """The method to draw the ball. Note: if you don’t
        specify the color of the ball, the color of each new
        element that is plotted will be randomly assigned.
        """
        #print(self.position)
        self.canvas.ax.plot(*self.position,"o", markersize=self.size)#, c=self.color) # moves with xpos = self.position and ypos = 0 (i.e. moves in 1D)

    def ZMF_collision_calculation(self, other):
            # a long calculation so I'll give it its own function
            # could probably use matrices

            # Initialisation
            u1, u2 = self.velocity, other.velocity
            m1, m2 = self.mass, other.mass

            # speed of ZMF frame:
            V_zmf = (m1*u1+m2*u2)/(m1+m2)
            
            # now define new velocities in ZMF
            u1_zmf, u2_zmf = u1 - V_zmf, u2 - V_zmf

            ### Now find the proportions of it in the r and r_perp directions (where r is the displacement between radii)
            # first find r and r_perp
            r = other.position - self.position
            rMag = np.linalg.norm(r)
            #...and normalising
            r = r / rMag

            r_perp = np.zeros(2)
            r_perp[0], r_perp[1] = -r[1], r[0]
            #print(f"dot = {np.dot(r,r_perp)}") # --> always 0 so they're always perpendicular!

            # then dot u_zmf with r and r_perp
            # a and b are constants
            a1, b1 = np.dot(u1_zmf, r), np.dot(u1_zmf, r_perp)
            a2, b2 = np.dot(u2_zmf, r), np.dot(u2_zmf, r_perp)

            # now we can find v_zmf's
            v1_zmf = b1 * r_perp - a1 * r
            v2_zmf = b2 * r_perp - a2 * r

            #and now we can find v1 and v2, out of the ZMF
            v1, v2 = v1_zmf + V_zmf, v2_zmf + V_zmf

            return v1, v2




    def collide(self, other):
        # Watch out with the threshold, this should
        # be greater than velocity. If, in a single
        # iteration, balls move relatively more than
        # what the threshold is, the collision might
        # not be triggered.

        #### we can define the threshold depending on the speed and radius of the particles
        sum_of_radii = self.radius + other.radius
        #speed_of_approach = abs(self.velocity-other.velocity)
        #so in one frame,    changes by -speed_of_appr * dT (if they are moving together)


          
        threshold = sum_of_radii# - speed_of_approach*self.canvas.dT
        #^ the effect of speed_of_approach on the threshold I am unsure of.
        if np.linalg.norm(self.position - other.position) < threshold: #just x direction
            # Need to use ZMF for this


            v1, v2 = self.ZMF_collision_calculation(other)
            # u1, u2 = self.velocity, other.velocity
            # m1, m2 = self.mass, other.mass
            # v1 = ((1-m2/m1)/(1+m2/m1))*u1 + (2/(1+m1/m2))*u2
            # v2 = (2/(1+m2/m1))*u1 - ((1-m2/m1)/(1+m2/m1))*u2

            
            self.velocity = v1
            other.velocity = v2

        # Wall collision
        for ball in [self, other]:
            collision = False
            #assert abs(ball.position[0]) < self.canvas.size/2 - ball.radius and abs(ball.position[1]) < self.canvas.size/2 - ball.radius
            # if abs(ball.position[0]) > self.canvas.size/2 - ball.radius:
            #     print("x-bounce")
            #     ball.velocity[0] *= -1 #bounce off the vertical wall if radius touches it

            if ball.position[0] > self.canvas.size/2 - ball.radius:
                #print("right x-bounce")
                old_velocity = ball.velocity[0]
                ball.velocity[0] = -abs(ball.velocity[0]) #bounce off the vertical wall if radius touches it
                new_velocity = ball.velocity[0]
                collision = True
            elif -1*ball.position[0] > self.canvas.size/2 - ball.radius:
                #print("left x-bounce")
                old_velocity = ball.velocity[0]
                ball.velocity[0] = abs(ball.velocity[0]) #bounce off the vertical wall if radius touches it
                new_velocity = ball.velocity[0]
                collision = True

            if ball.position[1] > self.canvas.size/2 - ball.radius:
                #print("top y-bounce")
                old_velocity = ball.velocity[1]
                ball.velocity[1] = -abs(ball.velocity[1]) #bounce off the horizontal wall if radius touches it           
                new_velocity = ball.velocity[1]
                collision = True
            elif -1*ball.position[1] > self.canvas.size/2 - ball.radius:
                #print("bottom y-bounce")
                old_velocity = ball.velocity[1]
                ball.velocity[1] = abs(ball.velocity[1]) #bounce off the horizontal wall if radius touches it
                new_velocity = ball.velocity[1]
                collision = True
            
            if collision:
                momentum_change = 2 * ball.mass * abs(old_velocity - new_velocity)
                self.canvas.pressure += momentum_change



#Main body of code

canvas = Canvas()

#Initialising all the balls
ball1 = Ball(canvas, mass=1, position=[-2,0], velocity=[0.05,0.01])
ball2 = Ball(canvas, mass=4, position=[2,0], velocity=[-0.07,-0.05])
ball3 = Ball(canvas, mass=5, position=[4,0], velocity=[-0.03,0.01])
ball4 = Ball(canvas, mass=2, position=[-5,0], velocity=[0.05,0.1])
ball5 = Ball(canvas, mass=3, position=[8,0], velocity=[-0.02,0.01])
ball6 = Ball(canvas, mass=3, position=[8,5], velocity=[-0.02,0.01])
ball7 = Ball(canvas, mass=4, position=[-4,4], velocity=[0.02,-0.04])
ball8 = Ball(canvas, mass=3, position=[7,-7], velocity=[-0.02,0.01])
ball9 = Ball(canvas, mass=1, position=[0,2], velocity=[0.05,0.1])
ball10 = Ball(canvas, mass=4, position=[2,5], velocity=[-0.07,-0.05])
ball11 = Ball(canvas, mass=5, position=[-3,8], velocity=[-0.03,0])
ball12 = Ball(canvas, mass=2, position=[-9,0], velocity=[0.05,0.1])
ball13 = Ball(canvas, mass=3, position=[-8,0], velocity=[-0.02,0])
ball14 = Ball(canvas, mass=3, position=[8,-5], velocity=[-0.02,0.01])
ball15 = Ball(canvas, mass=2, position=[-4,-4], velocity=[0.02,-0.04])
ball16 = Ball(canvas, mass=3, position=[7,8], velocity=[-0.02,0.01])
ball17 = Ball(canvas, mass=3, position=[9,9], velocity=[-0.02,-0.03])
ball18 = Ball(canvas, mass=3, position=[-7,4], velocity=[-0.05,0.01])
ball19 = Ball(canvas, mass=2, position=[-4,-6], velocity=[0.02,-0.06])
ball20 = Ball(canvas, mass=3, position=[9,7], velocity=[-0.02,0.01])

def animate(i):
    """This controls the animation."""
    print("The frame is:", i)
    canvas.update_balls()
    canvas.check_collision()
    canvas.fix_axes()
    if i % 10 == 1:
        canvas.total_kinetic_energy_calc()
        print(f"KE = {canvas.total_KE}")
        if i % 110 == 101:
            print("As you can see, the total KE of the system is constant!!")

        # average_pressure = canvas.pressure / (i+1)
        # print(canvas.pressure)
        ### ^ This was not very useful, so I've decided not to include it.


# This calls the animate function and creates animation.
anim = animation.FuncAnimation(canvas.fig, animate, frames=250, interval=10)
plt.show()

# # This saves the animation.
#anim.save("Session6_Collisions/balls2d_animation_ZMF.gif", writer="pillow", dpi=200)
