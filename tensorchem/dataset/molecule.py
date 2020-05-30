"""
Dataset objects for holding molecular geometries (atomic numbers, coordinates)
and their properties (energy, atomic forces, dipole moments).
"""

import json
import torch


class MoleculeSet:
    def __init__(self, atomic_nums=None):
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
            "atomic_numbers": self.atomic_nums,
            "geometries": [geom.export_json() for geom in self.geometries],
            "identifiers": self.identifiers
        }
        print(json_data)
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
        self.atomic_nums = json_data['atomic_numbers']
        self.identifiers = json_data['identifiers']
        for geom in json_data['geometries']:
            self.geometries.append(Geometry.from_json(self.atomic_nums, geom))


class Geometry:
    def __init__(self, atomic_nums=(0,), coords=((0.0, 0.0, 0.0),), properties={}):
        self.atomic_nums = atomic_nums
        self.coords = coords
        self.properties = properties

    def __repr__(self):
        return str(self.export_json())

    def __str__(self):
        rep = str(self.n_atoms)
        rep += "\n\n"
        for i, at_num in enumerate(self.atomic_nums):
            rep += "     ".join((str(at_num), str(self.coords[i][0]), str(self.coords[i][1]), str(self.coords[i][2])))
            rep += "\n"
        return rep

    @classmethod
    def from_json(cls, atomic_nums, json_dict):
        new_geom = cls()
        new_geom.atomic_nums = atomic_nums
        new_geom.coords = json_dict['coordinates']
        new_geom.properties = json_dict['properties']
        return new_geom

    def export_json(self):
        return {"atomic_nums": self.atomic_nums,
                "coordinates": self.coords,
                "properties": self.properties}

    @property
    def n_atoms(self):
        return len(self.atomic_nums)

    def __get_dist_matrix__(self):
        return torch.norm(self.coords)


if __name__ == "__main__":
    mol1 = MoleculeSet()
    mol1.load('19021.mset')
