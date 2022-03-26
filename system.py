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

    def run(self):
        for body in self.smallBodyList:
            body.updatePosition(self.largeBodyList, self.timeStep)

    def plotAllPositions(self):

