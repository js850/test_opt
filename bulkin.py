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
    target_E = 0
    def getEnergy(self, coords):
        x, y = coords
        E = 100. * sqrt(abs(y - 0.01 * x**2)) + 0.01 * abs(x + 10.)
        return E
    
    def get_limits(self):
        xmin = np.array([-15,-3])
        xmax = np.array([-5,3])
        return xmin, xmax



class BulkinSystem(BaseSystem):    
    def get_potential(self):
        return Bulkin()
    
    def get_random_configuration(self):
        x = np.random.uniform(-13,-8, 1)
        y = np.random.uniform(-4,-4, 1)
        return np.array([x,y]).flatten()
         

if __name__ == "__main__":
    from ackey import BenchmarkSystem
    mysys = BenchmarkSystem(Bulkin())
    mysys.do_benchmark()
    exit(1)
    

    system = BulkinSystem()
    print system.get_random_configuration()
    print system.get_random_minimized_configuration()
    
    db = system.create_database()
    bh = system.get_basinhopping(db)
    for i in range(20):
        bh.run(1)
        print bh.markovE, bh.coords
    
    print "minima found"
    for m in db.minima():
        print m.energy, m.coords
