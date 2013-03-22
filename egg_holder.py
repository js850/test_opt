from numpy import exp, cos, sqrt, pi, sin
import numpy as np

from base_function import BaseFunction

#class EHBenchmarkSystem(BenchmarkSystem):
#    def __init__(self, *args):
#        super(EHBenchmarkSystem, self).__init__(*args)
#        self.params.basinhopping.temperature = 10.
#    def get_takestep(self, **kwargs):
#        print "mytakestep is being called"
#        if "stepsize" in kwargs:
#            kwargs.pop("stepsize")
#        return super(EHBenchmarkSystem, self).get_takestep(stepsize=100., verbose=True, interval=50., **kwargs)
        

class EggHolder(BaseFunction):
    xmin = np.array([-512.,-512.])
    xmax = np.array([512.,512.])
    target_E = -959.6407
    def getEnergy(self, coords):
        x, y = coords
        E = -(y+47.) * sin(sqrt(abs(y + 0.5*x + 47.))) - x * sin(sqrt(abs(x - (y + 47.))))
        return E
    
    def dabs(self, x):
        if x > 0: return 1.
        elif x < 0: return -1.
        else: return 0.
    
    def getEnergyGradient_2(self, coords):
        """there is a bug here.  it doesn't pass the test"""
        E = self.getEnergy(coords)
        x, y = coords
        dEdy = (-sin(sqrt(abs(47. + 0.5*x + y))) + 
            (x*cos(sqrt(abs(-47. + x - y)))
             * -1.*self.dabs(-47. + x - y))/
            (2.*sqrt(abs(-47. + x - y))) 
            + 
            (0.5*(-47. - y)*cos(sqrt(abs(47. + 0.5*x + y)))
             * self.dabs(47. + 0.5*x + y)
            )/sqrt(abs(47. + 0.5*x + y)))
        
        dEdx = (
            -sin(sqrt(abs(-47. + x - y))) 
            - 
            (x*cos(sqrt(abs(-47. + x - y)))
             *self.dabs(-47. + x - y))/
            (2.*sqrt(abs(-47. + x - y))) 
            + 
            (0.25*(-47. - y)*cos(sqrt(abs(47. + 0.5*x + y)))
            *0.5*self.dabs(47. + 0.5*x + y)
            )/sqrt(abs(47. + 0.5*x + y)))
        grad = np.array([dEdx, dEdy])
        print "grad", grad
        return E, grad

if __name__ == "__main__":
    f = EggHolder()
#    f.test_potential(f.get_random_configuration())
    
    from base_function import makeplot2d
    makeplot2d(EggHolder())
    