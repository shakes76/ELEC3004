# -*- coding: utf-8 -*-
"""
Simulate the trajectory of a craft around the earth in 2D using Newton's laws of motion, 
but using Euler's method

https://astronomy.stackexchange.com/questions/7806/exercise-2d-orbital-mechanics-simulation-python

Most notably:
If you apply something simple like forward Euler or backward Euler, you will see 
the Earth spiral out to infinity or in toward the sun, respectively, but that is 
an effect of the numerical errors. If you use a more accurate method, like the 
classical 4th-order Runge-Kutta method, you'll find that it stays close to a true 
orbit for a while but still eventually goes off to infinity. The right approach is 
to use a symplectic method, which will keep the Earth in the correct orbit -- 
although its phase will still be off due to numerical errors.

Created on Tue Apr 24 14:11:40 2018

@author: uqscha22
"""
import matplotlib.pyplot as plt
import math
plt.ion()

G = 6.673e-11  # gravity constant
center = 100.0
gridArea = [0, 2*center, 0, 2*center]  # margins of the coordinate grid
gridScale = 10000000  # 1 unit of grid equals 10000000m or 10000km

fig, ax = plt.subplots(figsize=(8, 8))

plt.clf()  # clear plot area
plt.axis(gridArea)  # create new coordinate grid
plt.grid(b="on")  # place grid

class Object:
    _instances = []
    def __init__(self, name, position, radius, mass, color="black"):
        self.name = name
        self.position = position
        self.radius = radius  # in grid values
        self.mass = mass
        self.color = color
        self.placeObject()
        self.velocity = 0
        Object._instances.append(self)

    def placeObject(self):
        drawObject = plt.Circle(self.position, radius=self.radius, fill=False, color=self.color)
        plt.gca().add_patch(drawObject)
        plt.show()

    def giveMotion(self, deltaV, motionDirection, time):
        if self.velocity != 0:
            x_comp = math.sin(math.radians(self.motionDirection))*self.velocity
            y_comp = math.cos(math.radians(self.motionDirection))*self.velocity
            x_comp += math.sin(math.radians(motionDirection))*deltaV
            y_comp += math.cos(math.radians(motionDirection))*deltaV
            self.velocity = math.sqrt((x_comp**2)+(y_comp**2))

            if x_comp > 0 and y_comp > 0:  # calculate degrees depending on the coordinate quadrant
                self.motionDirection = math.degrees(math.asin(abs(x_comp)/self.velocity))  # update motion direction
            elif x_comp > 0 and y_comp < 0:
                self.motionDirection = math.degrees(math.asin(abs(y_comp)/self.velocity)) + 90
            elif x_comp < 0 and y_comp < 0:
                self.motionDirection = math.degrees(math.asin(abs(x_comp)/self.velocity)) + 180
            else:
                self.motionDirection = math.degrees(math.asin(abs(y_comp)/self.velocity)) + 270

        else:
            self.velocity = self.velocity + deltaV  # in m/s
            self.motionDirection = motionDirection  # degrees
        self.time = time  # in seconds
        self.vectorUpdate()

    def vectorUpdate(self):
        self.placeObject()
        data = []

        for t in range(self.time):
            motionForce = self.mass * self.velocity  # F = m * v
            x_net = 0
            y_net = 0
            for x in [y for y in Object._instances if y is not self]:
                distance = math.sqrt(((self.position[0]-x.position[0])**2) +
                             (self.position[1]-x.position[1])**2)
                gravityForce = G*(self.mass * x.mass)/((distance*gridScale)**2)

                x_pos = self.position[0] - x.position[0]
                y_pos = self.position[1] - x.position[1]

                if x_pos <= 0 and y_pos > 0:  # calculate degrees depending on the coordinate quadrant
                    gravityDirection = math.degrees(math.asin(abs(y_pos)/distance))+90

                elif x_pos > 0 and y_pos >= 0:
                    gravityDirection = math.degrees(math.asin(abs(x_pos)/distance))+180

                elif x_pos >= 0 and y_pos < 0:
                    gravityDirection = math.degrees(math.asin(abs(y_pos)/distance))+270

                else:
                    gravityDirection = math.degrees(math.asin(abs(x_pos)/distance))

                x_gF = gravityForce * math.sin(math.radians(gravityDirection))  # x component of vector
                y_gF = gravityForce * math.cos(math.radians(gravityDirection))  # y component of vector

                x_net += x_gF
                y_net += y_gF

            x_mF = motionForce * math.sin(math.radians(self.motionDirection))
            y_mF = motionForce * math.cos(math.radians(self.motionDirection))
            x_net += x_mF
            y_net += y_mF
            netForce = math.sqrt((x_net**2)+(y_net**2))

            if x_net > 0 and y_net > 0:  # calculate degrees depending on the coordinate quadrant
                self.motionDirection = math.degrees(math.asin(abs(x_net)/netForce))  # update motion direction
            elif x_net > 0 and y_net < 0:
                self.motionDirection = math.degrees(math.asin(abs(y_net)/netForce)) + 90
            elif x_net < 0 and y_net < 0:
                self.motionDirection = math.degrees(math.asin(abs(x_net)/netForce)) + 180
            else:
                self.motionDirection = math.degrees(math.asin(abs(y_net)/netForce)) + 270

            self.velocity = netForce/self.mass  # update velocity
            traveled = self.velocity/gridScale  # grid distance traveled per 1 sec
            self.position = (self.position[0] + math.sin(math.radians(self.motionDirection))*traveled,
                             self.position[1] + math.cos(math.radians(self.motionDirection))*traveled)  # update pos
            data.append([self.position[0], self.position[1]])

            collision = 0
            for x in [y for y in Object._instances if y is not self]:
                if (self.position[0] - x.position[0])**2 + (self.position[1] - x.position[1])**2 <= x.radius**2:
                    collision = 1
                    break
            if collision != 0:
                print("Collision!")
                break

        plt.plot([x[0] for x in data], [x[1] for x in data])

Earth = Object(name="Earth", position=(center, center), radius=6.371, mass=5.972e24, color='black')

distance = 360000*1000/gridScale
Moon = Object(name="Moon", position=(center+distance, center+distance), radius=1.737, mass = 7.347e22, color='green')  # position not to real scale

totalTime = 29*24*60*60 #one lunar month
speed = 900 #1.022 km/s
Moon.giveMotion(deltaV=speed, motionDirection=100, time=int(4*totalTime))

plt.show(block=True)