import numpy as np
from numpy import exp, sqrt, cos, pi, sin

from base_function import BaseFunction

class Ackley(BaseFunction):
    #note: this function is not smooth at the origin.  the gradient will never
    #converge in the minimizer
    target_E = 0.
    xmin = np.array([-5.,-5.])
    xmax = np.array([5.,5.])
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
    
if __name__ == "__main__":
    f = Ackley()
    f.test_potential(f.get_random_configuration())
    
    from base_function import makeplot2d
    makeplot2d(f)

