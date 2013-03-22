"""
http://en.wikipedia.org/wiki/Test_functions_for_optimization

minimum: f(-10,1) = 0 
for 
-15 < x < -5
-3 < y < 3
"""


from numpy import exp, cos, sqrt, pi
import numpy as np

from base_function import BaseFunction
from pygmin.systems import BaseSystem

class Bulkin(BaseFunction):
    xmin = np.array([-15.,-3.])
    xmax = np.array([-5.,3.])
    target_E = 0.
    def getEnergy(self, coords):
        x, y = coords
        E = 100. * sqrt(abs(y - 0.01 * x**2)) + 0.01 * abs(x + 10.)
        return E

    def get_random_configuration(self):
        x = np.random.uniform(-13,-8, 1)
        y = np.random.uniform(-4,-4, 1)
        return np.array([x,y]).flatten()


class BukinSystem(BaseSystem):    
    def get_potential(self):
        return Bulkin()
    
         

if __name__ == "__main__":
    from base_function import makeplot2d
    makeplot2d(Bulkin())
