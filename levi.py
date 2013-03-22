import numpy as np
from numpy import exp, sqrt, cos, pi, sin

from base_function import BaseFunction


class Levi(BaseFunction):
    target_E = 0.
    xmin = np.array([-10.,-10.])
    xmax = np.array([10.,10.])

    def getEnergy(self, coords):
        x, y = coords
        E = sin(3.*pi*x)**2 + (x-1.)**2 * (1. + sin(3*pi*y)**2) \
            + (y-1.)**2 * (1. + sin(2*pi*y)**2)
        return E
    
    def getEnergyGradient(self, coords):
        x, y = coords
        E = self.getEnergy(coords)
        
        dEdx = 2.*3.*pi* cos(3.*pi*x) * sin(3.*pi*x) + 2.*(x-1.) * (1. + sin(3*pi*y)**2)
        
        dEdy = (x-1.)**2 * 2.*3.*pi* cos(3.*pi*y) * sin(3.*pi*y) + 2. *  (y-1.) * (1. + sin(2*pi*y)**2) \
            + (y-1.)**2 * 2.*2.*pi * cos(2.*pi*y) * sin(2.*pi*y)
        
        return E, np.array([dEdx, dEdy])

if __name__ == "__main__":
    f = Levi()
    f.test_potential(f.get_random_configuration())
    
    from base_function import makeplot2d
    makeplot2d(f)

