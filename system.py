"""
The system class allows for the management and observation of a set of bodies
making up a solar system.

Eventually one will be able to observe the trajectory of a small body moving
through this system, as affected by the large bodies.
"""


import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

import body


class System():
    def __init__(self, largeBodyFile):

        #create the large bodies from the text file
        self.largeBodyList = []
        #get all the planet data as strings
        largeBodyData = np.loadtxt(largeBodyFile, skiprows=1)

        #gather a list of the large body objects
        for row in range(len(largeBodyData)):
            currentName = largeBodyData[row, 0]
            currentParent = largeBodyData[row, 1]
            currentMass = largeBodyData[row, 2]
            currentRadius = largeBodyData[row, 3]
            current_e = largeBodyData[row, 4]
            currentSA = largeBodyData[row, 5]
            currentInc = largeBodyData[row, 6]
            currentLAN = largeBodyData[row, 7]
            currentAP = largeBodyData[row, 8]
            currentMA = largeBodyData[row, 9]

            #TODO fix the mess with setting things to moons properly etc

            if currentParent != "na":
                #need to select the parent body, which should already exist
                #thanks to the order in the bodyData file.
                found = False
                counter = 0
                while counter < len(self.largeBodyList):
                    if self.largeBodyList[counter].name == currentParent:
                        currentParent = self.largeBodyList[counter]
                        found = True

                if not found:
                    raise(Exception("Parent body of " + currentName + " not found!"))

     
            self.largeBodyList.append(body.largeBody(row, currentName, currentParent, currentMass, currentRadius, current_e, currentSA, currentInc, currentLAN, currentAP, currentMA))

        # self.smallBodyList = smallBodyList

        # self.positionHistories = []

        self.animationTimeStep = 60*60*24 #a day's worth of seconds
        return

    def run(self, N, timeStep):
        """
        N:          The number of steps to run for, integer.
        timeStep:   The amount of time between each step, seconds, float.
        """
        counter = 0

        while counter < N:

            for body in self.largeBodyList:
                self.positionHistories[body.ID].append(body.updatePosition(timeStep))
                
            # for body in self.smallBodyList:
            #     Íself.positionHistories[body.ID].append(body.updatePosition(self.largeBodyList, timeStep)[0])

            counter += 1

        self.positionHistories = np.array(self.positionHistories)
        
        return

    def getLargeBodyPositions(self, absoluteTime):
        """
        Get the positions of the large bodies at a particular time.
        Return as an array in the same order as the planets in order
        of ID.
        """
        numLargeBodies = len(self.largeBodyList)
        positionsList = np.zeros((numLargeBodies, 3))
        for i in range(0, numLargeBodies):
            positionsList[i] = body.getPosition(absoluteTime)




        return

    def animateSystem(i):
        return



    def plotPositions(self, positionsList, figureObject):

        
        for i in range(0, np.shape(self.positionHistories)[0]):
            plt.plot(self.positionHistories[i, :, 0], self.positionHistories[i, :, 1])
        
        return

