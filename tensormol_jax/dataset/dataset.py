"""
Dataset objects for holding molecular geometries (atomic numbers, coordinates)
and their properties (energy, atomic forces, dipole moments) for training models.
"""

import json
import numpy as np 

from torch.utils.data import Dataset as TorchDataset
from tensormol_jax.dataset.molecule import MoleculeSet, Geometry

class Dataset(TorchDataset):
    def __init__(self):
        super(Dataset, self).__init__()
        self.unique_atoms = []
        self.filename = None

    def _init_dataset(self):
        return


class MolDataset(Dataset):
    def __init__(self):
        super(MolDataset, self).__init__()
        self.molecules = () # Immutable type so order of molecules cannot change during training
        self.idx_map = {} # Maps an overall sample index to the molecule and geometry indices

    def __len__(self):
        return sum([len(molecule_set) for molecule_set in self.molecules])

    def __getitem__(self, idx):
        return


class MixedDataset(Dataset):
    def __init__(self):
        super(MixedDataset, self).__init__()
        self.samples = []

    def __len__(self):
        return sum([len(molecule_set) for molecule_set in self.molecules])

    def __getitem__(self, idx):
        return

    def save(self, filename=None):
        if filename is None:
            if self.filename is None:
                print("No filename given for saving")
                exit(0)
            else:
                filename = self.filename
        json_data = {
            "atomic_number": [sample['atomic_number'] for sample in self.samples],
            "coordinates": [geom.coords.tolist() for sample in self.samples for geom in sample['geometries']],
            "properties": [geom.properties for sample in self.samples for geom in sample['geometries']],
        }
        with open(filename, "w") as f:
            json.dump(json_data, f)

    def load(self, filename=None):
        if filename is None:
            try:
                filename = self.filename
            except:
                print("No file specified for Dataset loading.")
        with open(filename, 'r') as f:
            for line in f:
                mset = json.loads(line)
                sample = {
                    "atomic_number": mset['atomic_number'],
                    "geometries": [Geometry(np.array(mset['coordinates']), mset['properties'])]
                }
                self.samples.append(sample)


if __name__ == "__main__":
    ds1 = MixedDataset()
    ds1.load('/home/adriscoll/tensormol-jax/tensormol_jax/data/ani1x-mol.mset')
    ds1.save('/home/adriscoll/tensormol-jax/tensormol_jax/data/ani1x-saved.mset')
