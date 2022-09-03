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
    
    def getPosition(self, format="normal"):
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
    def __init__(self, ID, bodyName, bodyParent, mass, radius, eccentricity, semimajorAxis, inclination, longAscNode, argPeriapsis, meanAnomaly):
        """
        A largeBody moves on 'rails'.  For this, we do not need the state vector,
        just the Keplerian orbital elements.
        """
        # super().__init__(ID, mass, position, velocity)
        #body parameters
        self.ID = ID
        self.name = bodyName
        self.parent = bodyParent
        self.m = mass
        self.r = radius
        #orbital parameters
        self.e = eccentricity
        self.a = semimajorAxis
        self.inc = np.deg2rad(inclination) #these two are stored as degrees
        self.longAN = np.deg2rad(longAscNode) #because thats what ksp does in the wiki
        self.argPer = argPeriapsis
        self.meanAnom = meanAnomaly

        #Compute some useful values for later.  See "Position as a function of time"
        #on wikipedia: https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Position_as_a_function_of_time
        self.T = 0 #compute the period somehow
        
        self.n = 2*np.pi/T


        return 

    
    def keplersEquation(E, M, epsilon):
        #This is Kepler's equation (rearranged)
        # See https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Position_as_a_function_of_time
        return M - epsilon*np.sin(E)

    def keplersDerivative(E, M, epsilon):
        #derivative of Kepler's equation
        return -epsilon*np.cos(E)

    

    def getPosition(self, t):
        #t is the time since the beginning of the simulation
        currentMeanAnom = self.n*t + self.meanAnom

        #calculate eccentric anomaly E
        E = sp.optimize.newton(self.keplersEquation, np.pi/4, self.keplersDerivative, args=(currentMeanAnom, self.e))
        #calculate the true anomaly theta
        #TODO calculate this
        #calculate the heliocentric distance r
        r = self.a*(1-self.e*np.cos(E))


        return self.r, self.v

