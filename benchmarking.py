"""
benchmarking basinhopping for functions from 

http://en.wikipedia.org/wiki/Test_functions_for_optimization
"""
from numpy import exp, cos, sqrt, pi, sin
import numpy as np

from base_function import BaseFunction
from pygmin.systems import BaseSystem

class BenchmarkSystem(BaseSystem):
    Etol = 1e-6
    def __init__(self, potential):
        BaseSystem.__init__(self)
        self.pot = potential
        try:
            self.params.basinhopping.temperature = self.pot.temperature
        except AttributeError:
            pass
    
    def test_potential(self):
        coords = self.get_random_configuration()
        e = self.pot.getEnergy(coords)
        egrad, grad = self.pot.getEnergyGradient(coords)
        gradnum = self.pot.NumericalDerivative(coords, 1e-6)
        print "testing energy"
        print e
        print egrad
        print np.abs(e-egrad)

        print "testing gradient"
        print grad
        print gradnum
        print np.max(np.abs(grad-gradnum))
    
    def get_potential(self):
        return self.pot
    
    def get_random_configuration(self):
        if hasattr(self.pot, "get_random_configuration"):
            return self.pot.get_random_configuration()
        xmin, xmax = self.pot.xmin, self.pot.xmax
        x = np.random.uniform(xmin[0] + .01, xmax[0] - .01)
        y = np.random.uniform(xmin[1] + .01, xmax[1] - .01)
        return np.array([x,y])
    
    def accept_test(self, E, coords, *args, **kwargs):
        if not hasattr(self.pot, "xmin"): return True
        if np.any(coords < self.pot.xmin):
            return False
        if np.any(coords > self.pot.xmax):
            return False
        return True
    
    def stop_criterion(self, coords, E, *args, **kwargs):
        if E < self.pot.target_E + self.Etol:
            return True
        else:
            return False
        
    def get_takestep(self):
        if hasattr(self.pot, "stepsize"):
            return BaseSystem.get_takestep(self, stepsize=self.pot.stepsize)
        else:
            return BaseSystem.get_takestep(self)
    
    def do_benchmark(self):
        x0 = self.get_random_configuration()
        bh = self.get_basinhopping(confCheck = [self.accept_test], outstream=None)
        for i in range(1000):
            bh.run(1)
            if self.stop_criterion(bh.coords, bh.markovE):
                print "done"
                print "found the global minimum after", i, "basinhopping steps"
                print "E", bh.markovE
                print "coords", bh.coords
                break
    
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


from base_function import LJ
class LJ30(LJ):
    target_E = -128.286571
    natoms = 30
    def get_random_configuration(self):
        return np.random.uniform(-1,1,[3*self.natoms]) * float(self.natoms)**(1./3) 

class Schaffer2(BaseFunction):
    target_E = 0.
    xmin = np.array([-100,-100])
    xmax = np.array([100,100])


if __name__ == "__main__":
    print ""
    print "doing benchmark for Ackey function"
    mysys = BenchmarkSystem(Ackey())
    mysys.test_potential()
#    exit(1)
    mysys.do_benchmark()

    print ""
    print "doing benchmark for Levi function"
    mysys = BenchmarkSystem(Levi())
    mysys.test_potential()
#    exit(1)
    mysys.do_benchmark()

    print ""
    print "doing benchmark for Holder Table function"
    mysys = BenchmarkSystem(HolderTable())
    mysys.test_potential()
##    exit(1)
    mysys.do_benchmark()

    print ""
    print "doing benchmark for LJ30 function"
    mysys = BenchmarkSystem(LJ30())
#    mysys.test_potential()
##    exit(1)
    mysys.do_benchmark()
