"""
Dataset objects for holding molecular geometries (atomic numbers, coordinates)
and their properties (energy, atomic forces, dipole moments).
"""

import json
import torch
import numpy as np


class MoleculeSet:
    def __init__(self):
        super(MoleculeSet, self).__init__()
        self.atomic_nums = ()
        self.geometries = []
        self.identifiers = {}
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
                raise FileNotFoundError("No filename given for saving")
            else:
                filename = self.filename
        json_data = {
            "atomic_number": self.atomic_nums,
            "coordinates": [geom.coords.tolist() for geom in self.geometries],
            "properties": {key: value for prop in self.geometries for key, value in prop.properties.items()},
            "identifiers": {key: value for (key, value) in self.identifiers.items()}
        }
        with open(filename, 'w') as f:
            json.dump(json_data, f)

    def load(self, filename=None):
        if filename is None:
            if self.filename is None:
                raise FileNotFoundError("No filename given for loading")
            else:
                filename = self.filename
        with open(filename) as f:
            json_data = json.load(f)
        for key in json_data.keys():
            if type(json_data[key]) == dict:
                self.identifiers = {key: json_data[key]}
        self.atomic_nums = json_data['atomic_number']
        self.identifiers = json_data['identifiers']
        for i, coords in enumerate(json_data['coordinates']):
            properties = {key: value[i] for key, value in json_data['properties'].items()}
            self.geometries.append(Geometry(self.atomic_nums, coords, properties))


class Geometry:
    def __init__(self, atomic_nums, coords, properties=None):
        self.atomic_nums = ()
        self.coords = coords
        if properties is None:
            self.properties = {}
        else:
            self.properties = properties

    def __repr__(self):
        rep = str(self.n_atoms)
        rep += "\n\n"
        for i, at_num in enumerate(self.atomic_nums):
            rep += "     ".join((str(at_num), str(self.coords[i, 0]), str(self.coords[i, 1]), str(self.coords[i, 2])))
            rep += "\n"
        return rep

    @property
    def n_atoms(self):
        return len(self.atomic_nums)

    def __get_dist_matrix__(self):
        return torch.norm(self.coords)


if __name__ == "__main__":
    mol1 = MoleculeSet()
    mol1.load('../data/ani1x-mol.mset')
