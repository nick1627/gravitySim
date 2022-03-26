import numpy as np
import scipy as sp




class Body():
    def __init__(self, ID, mass, position, velocity):
        """
        A general astronomical body

        mass:       The mass of the object (kg)
        position:   The position of the object (m)
        velocity:   The velocity of the body (m/s)
        ID:         An identification number.  Should be integer, unique.
        """
        self.m = mass
        self.r = position
        self.v = velocity
        self.ID = ID

        return

    def updatePosition(self, largeBodyList, timeStep):
        """
        largeBodyList:      A list of largeBody objects
        """

        totalAcceleration = 0

        for body in largeBodyList:
            if body.ID != self.ID:
                relativePosition = body.r - self.r
                totalAcceleration += (body.m*relativePosition)/(np.linalg.norm(relativePosition))**3
        
        totalAccelration *= sp.constants.G

        self.v += totalAcceleration*timeStep
        self.r += self.v*timeStep

        return self.r, self.v
    
    def getPosition(self):
        return self.r

    def getVelocity(self):
        return self.v



class smallBody(Body):
    def __init__(self, ID, mass, position, velocity):
        """
        A small body moves under the field of the largeBodys
        """
        super().__init__(ID, mass, position, velocity)
        return

class largeBody(Body):
    def __init__(ID, mass, position, velocity):
        """
        A largeBody moves on 'rails'.
        """
        super().__init__(ID, mass, position, velocity)
        return 

