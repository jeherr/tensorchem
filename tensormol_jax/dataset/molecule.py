"""
Dataset objects for holding molecular geometries (atomic numbers, coordinates)
and their properties (energy, atomic forces, dipole moments).
"""

import json
import torch


class MoleculeSet:
    def __init__(self):
        super(MoleculeSet, self).__init__()
        self.atomic_nums = ()
        self.geometries = []
        self.identifiers = {}
        self.min_geom = None

    def __len__(self):
        return len(self.geometries)

    def __setitem__(self, idx, value):
        self.geometries[idx] = value

    def __getitem__(self, idx):
        return self.geometries[idx]

    def save(self, filename=None):
        # Will do some stuff to collect all necessary data into a format for JSON
        json_data = None
        with open(filename, "w") as f:
            json.dump(json_data, f)


class Geometry:
    def __init__(self, coords=None, properties=None):
        self.coords = coords
        self.properties = properties

    def get_dist_matrix(self):
        return torch.norm(self.coords)
