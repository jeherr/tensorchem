"""
Dataset objects for holding molecular geometries (atomic numbers, coordinates)
Dataset objects for holding molecular geometries (atomic numbers, coordinates)
and their properties (energy, atomic forces, dipole moments).
"""

import json
import torch
import numpy as np

class MoleculeSet():
    def __init__(self):
        super(MoleculeSet, self).__init__()
        self.atomic_nums = ()
        self.geometries = []
        self.identifiers = {}
        self.properties = {}
        self.min_geom = None
        self.filename = None

    def __len__(self):
        return len(self.geometries)

    def __setitem__(self, idx, value):
        self.geometries[idx] = Geometry(value)

    def __getitem__(self, idx):
        return self.geometries[idx]

    def save(self, filename=None):
        if filename is None:
            if self.filename is None:
                # TODO This should be implemented as an error
                print("No filename given for saving")
                exit(0)
            else:
                filename = self.filename
        print(type(self.atomic_nums))
        json_data = {
            "atomic_number": self.atomic_nums,
            "coordinates": [geom.coords.tolist() for geom in self.geometries],
            "properties": [geom.properties for geom in self.geometries],
            "identifiers": self.identifiers
        }
        with open(filename, 'w') as f:
            json.dump(json_data, f)

    def load(self, filename=None):
        if filename is None:
            try:
                filename = self.filename
            except:
                "No file specified for Molecule Set loading."
        with open(filename) as f:
            json_data = json.load(f)[0]
        self.atomic_nums = json_data['atomic_nums']
        self.geometries = [Geometry(geom) for geom in json_data['geometries']]
        # self.identifiers = {"identifiers": json_data['identifiers']}


class Geometry:
    def __init__(self, coords=None, properties=None):
        self.coords = coords
        self.properties = properties

    def __get_dist_matrix__(self):
        return torch.norm(self.coords)


if __name__ == "__main__":
    mol1 = MoleculeSet()
    mol1.load('../data/ani1x-release.json')
    mol1.save('../data/mol_test_save.json')
