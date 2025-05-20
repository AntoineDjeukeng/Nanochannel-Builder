import numpy as np
from math import sqrt

class UnitCell:
    def __init__(self, coords, d, cz):
        self.atoms = coords
        self.atom_types = []
        self.satom_types = []
        self.atom_names = []
        self.charges = []
        self.name = ''
        self.box = [d * sqrt(3), d * 3, cz]
        self.center = np.mean(coords, axis=0)

    def build(self, vectors):
        if vectors.shape[0] != 3:
            raise ValueError('vectors should be a 3-element vector')
        vec = vectors - self.center
        return self.atoms + vec
