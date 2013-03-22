"""
benchmarking Lennard-Jones clusters
"""
from base_function import BaseFunction
from base_function import LJ

class LJ30(LJ):
    target_E = -128.286571
    natoms = 30
    def get_random_configuration(self):
        return np.random.uniform(-1,1,[3*self.natoms]) * float(self.natoms)**(1./3) 
