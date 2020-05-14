"""
Dataset objects for holding molecular geometries (atomic numbers, coordinates)
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

    def save(self, filename='mol_data.json'):
	# Will do some stuff to collect all necessary data into a format for JSON
	json_data = None
	coordinates = []
	properties = []
	for geom in self.geometries:
		geometry = Geometry(geom)	
		coordinates.append(geometry.coords)
		properties.append(geometry.properties)
	json_data = {
        	"atomic_number":self.atomic_nums,
        	"coordinates":coordinates,
        	"properties":properties,
        	"identifiers":self.identifiers
        }
	with open(filename, 'w') as f:
		json.dump(json_data,f)


class Geometry:
    def __init__(self, coords=None, properties=None):
        self.coords = coords
        self.properties = properties

    def get_dist_matrix(self):
        return torch.norm(self.coords)

if __name__ == "__main__":
   with open('/home/adriscoll/tensormol-jax/tensormol_jax/data/ani1x-release.json') as x:
	json_mol_data = json.load(x)

   json_mol = json_mol_data[0]
   mol1 = MoleculeSet()
   mol1.atomic_nums = json_mol['atomic_numbers']
   mol1.geometries = json_mol['atomic_positions']
   mol1.save()
