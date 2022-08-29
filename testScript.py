import numpy as np
import scipy as sp
import matplotlib as plt

from system import System
from body import largeBody, smallBody

Earth = largeBody(0, 10**24, np.array([0, 0, 0]), np.array([0, 0, 0]))
Moon = smallBody(1, 0.0123*10**24, np.array([0.3844*10**9, 0, 0]), np.array([0, 1.022*10**3, 0]))

s1 = System([Earth], [Moon], 1000)
s1.run(100000)

s1.plotAllPositions()

plt.show()
