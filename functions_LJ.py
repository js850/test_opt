"""
benchmarking Lennard-Jones clusters
"""
import numpy as np

from base_function import BaseFunction
from pygmin.potentials import LJ

class LennardJones(LJ):
    """
    The Lennard Jones potential
    
    a mathematically simple model that approximates the interaction between a 
    pair of neutral atoms or molecules.    
    http://en.wikipedia.org/wiki/Lennard-Jones_potential
    
    E = sum_ij V(r_ij)
    
    where r_ij is the cartesian distance between atom i and atom j, and the
    pair potential has the form
    
    V(r) = 4 * eps * ( (sigma / r)**12 - (sigma / r)**6
    
    Notes
    -----
    the double loop over many atoms makes this *very* slow in Python.  If it
    were in a compiled language it would be much faster.
    """
    def get_random_configuration(self):
        return np.random.uniform(-1,1,[3*self.natoms]) * float(self.natoms)**(1./3) 


class LJ38(LennardJones):
    natoms = 38
    target_E = -173.928427

class LJ30(LennardJones):
    natoms = 30
    target_E = -128.286571

class LJ20(LennardJones):
    natoms = 20
    target_E = -77.177043

class LJ13(LennardJones):
    natoms = 13
    target_E = -44.326801
