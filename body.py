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
        self.i = np.deg2rad(inclination) #these two are stored as degrees
        self.longAN = np.deg2rad(longAscNode) #because thats what ksp does in the wiki
        self.argPer = argPeriapsis
        self.meanAnom = meanAnomaly

        if not self.parent == "na":
            #Compute some useful values for later.  See "Position as a function of time"
            #on wikipedia: https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Position_as_a_function_of_time
            self.GM = sp.constants.G*self.bodyParent.m
            #compute the period using Kepler's third law
            self.T = 2*np.pi*np.sqrt(self.a**3/self.GM)
            #now compute mean motion (mean angle per unit time)
            self.n = 2*np.pi/self.T

            #Here is a good place to calculate the transformation matrices too
            longANMatrix = np.array([[np.cos(self.longAN), -np.sin(self.longAN), 0], [np.sin(self.longAN), np.cos(self.longAN), 0], [0, 0, 1]])
            inclinationMatrix = np.array([[1, 0, 0], [0, np.cos(self.i), -np.sin(self.i)], [0, np.sin(self.i), np.cos(self.i)]])
            argPeriapsisMatrix = np.array([[np.cos(self.argPer), -np.sin(self.argPer), 0], [np.sin(self.argPer), np.cos(self.argPer), 0], [0, 0, 1]])
            
            self.transfMatrix = np.matmul(longANMatrix, np.matmul(inclinationMatrix, argPeriapsisMatrix))
            self.invTransfMatrix = np.linalg.inv(self.transfMatrix)

        return 

    
    def keplersEquation(E, M, epsilon):
        #This is Kepler's equation (rearranged)
        # See https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Position_as_a_function_of_time
        return M - epsilon*np.sin(E)

    def keplersDerivative(E, M, epsilon):
        #derivative of Kepler's equation
        return -epsilon*np.cos(E)

    

    def getPosition(self, t):
        if self.bodyParent != "na":
            #t is the time since the beginning of the simulation
            currentMeanAnom = self.n*t + self.meanAnom
            #calculate eccentric anomaly E
            E = sp.optimize.newton(self.keplersEquation, np.pi/4, self.keplersDerivative, args=(currentMeanAnom, self.e))
            #calculate the true anomaly theta
            theta = 2*np.sqrt((1+self.e)/(1-self.e))*np.tan(E/2)
            #calculate the heliocentric distance r
            r = self.a*(1-self.e*np.cos(E))
            #Now we have r and theta.  Need to convert this to a cartesian vector,
            #then convert that into the desired frame using the rest of the 
            #orbital elements.
            position = np.array([r*np.cos(theta), r*np.sin(theta), 0])
            #use inverse transformation matrix to convert to inertial frame
            position = np.matmul(self.invTransfMatrix, position)

        else:
            position = np.array([0, 0, 0])

        return position

