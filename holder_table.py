import numpy as np
from numpy import exp, sqrt, cos, pi, sin

from base_function import BaseFunction


class HolderTable(BaseFunction):
    target_E = -19.2085
    xmin = np.array([-10.,-10.])
    xmax = np.array([10.,10.])
    stepsize = 2.
    temperature = 2.
    def getEnergy(self, coords):
        x, y = coords
        E = - abs(sin(x)* cos(y) * exp(abs(1. - sqrt(x**2 + y**2)/ pi)))
        return E
    
    def dabs(self, x):
        """derivative of absolute value"""
        if x < 0: return -1.
        elif x > 0: return 1.
        else: return 0.

    def getEnergyGradient(self, coords):
        x, y = coords
        R = sqrt(x**2 + y**2)
        g = 1. - R / pi
        f = sin(x)* cos(y) * exp(abs(g))
        E = -abs(f)
        
        
        dRdx = x / R
        dgdx = - dRdx / pi
        dfdx = cos(x) * cos(y) * exp(abs(g)) + f * self.dabs(g) * dgdx
        dEdx = - self.dabs(f) * dfdx
        
        dRdy = y / R
        dgdy = - dRdy / pi
        dfdy = -sin(x) * sin(y) * exp(abs(g)) + f * self.dabs(g) * dgdy        
        dEdy = - self.dabs(f) * dfdy
        return E, np.array([dEdx, dEdy])

if __name__ == "__main__":
    f = HolderTable()
    f.test_potential(f.get_random_configuration())
    
    from base_function import makeplot2d
    makeplot2d(f)
