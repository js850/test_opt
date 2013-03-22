"""
benchmarking basinhopping for functions from 

http://en.wikipedia.org/wiki/Test_functions_for_optimization
"""
from numpy import exp, cos, sqrt, pi, sin
import numpy as np


class Ackey(BaseFunction):
    target_E = 0.
    xmin = np.array([-5,-5])
    xmax = np.array([5,5])
    def getEnergy(self, coords):
        x, y = coords
        E = -20. * exp(-0.2 * sqrt(x**2 + y**2)) - \
            exp(0.5 * (cos(2. *pi * x) + cos(2. * pi * y))) + 20. + np.e
        return E
    
    def getEnergyGradient(self, coords):
        E = self.getEnergy(coords)
        x, y = coords
        R = sqrt(x**2 + y**2)
        term1 = -20. * exp(-0.2 * R)
        term2 = -exp(0.5 * (cos(2. *pi * x) + cos(2. * pi * y)))
        
        deriv1 = term1 * (-0.2 * 0.5 / R)
        
        dEdx = 2.* deriv1 * x  - term2 * pi * sin(2.*pi*x)
        dEdy = 2.* deriv1 * y  - term2 * pi * sin(2.*pi*y)
        
        return E, np.array([dEdx, dEdy])

class Levi(BaseFunction):
    target_E = 0.
    xmin = np.array([-10,-10])
    xmax = np.array([10,10])
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
  
class HolderTable(BaseFunction):
    target_E = -19.2085
    xmin = np.array([-10,-10])
    xmax = np.array([10,10])
    stepsize = 2.
    temperature = 2.
    def getEnergy(self, coords):
        x, y = coords
        E = - abs(sin(x)* cos(y) * exp(abs(1. - sqrt(x**2 + y**2)/ pi)))
        return E
    
    def dabs(self, x):
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


class Schaffer2(BaseFunction):
    target_E = 0.
    xmin = np.array([-100,-100])
    xmax = np.array([100,100])
