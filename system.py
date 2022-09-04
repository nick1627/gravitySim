import numpy as np
import scipy as sp
import matplotlib as plt

import body


class System():
    def __init__(self, largeBodyFile, smallBodyList, timeStep):

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
     
            self.largeBodyList.append(body.largeBody(row, currentName, currentParent, currentMass, currentRadius, current_e, currentSA, currentInc, currentLAN, currentAP, currentMA))

        self.smallBodyList = smallBodyList
        self.timeStep = timeStep

        self.positionHistories = []
        return

    def run(self, N):
        """
        N:      The number of steps to run for
        """
        counter = 0

        while counter < N:

            for body in self.largeBodyList:
                self.positionHistories[body.ID].append(body.updatePosition(self.timeStep))
                
            for body in self.smallBodyList:
                self.positionHistories[body.ID].append(body.updatePosition(self.largeBodyList, self.timeStep)[0])

            counter += 1

        self.positionHistories = np.array(self.positionHistories)
        
        return

    def plotAllPositions(self):

        plt.figure(1)
        for i in range(0, np.shape(self.positionHistories)[0]):
            plt.plot(self.positionHistories[i, :, 0], self.positionHistories[i, :, 1])
        
        return

