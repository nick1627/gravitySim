import numpy as np
import scipy as sp
import matplotlib as plt

import body


class System():
    def __init__(self, largeBodyList, smallBodyList, timeStep):
        self.largeBodyList = largeBodyList
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

