"""
Dataset objects for holding molecular geometries (atomic numbers, coordinates)
and their properties (energy, atomic forces, dipole moments) for training models.
"""

import json
from torch.utils.data import Dataset as TorchDataset


class Dataset(TorchDataset):
    def __init__(self):
        super(Dataset, self).__init__()
        self.unique_atoms = []

    def _init_dataset(self):
        return

    def save(self, filename=None):
        # Will do some stuff to collect all necessary data into a format for JSON
        json_data = {
	   "unique_atoms":self.unique_atoms
	}
        with open(filename, "w") as f:
            json.dump(json_data, f)


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
        self.geometries = []

    def __len__(self):
        return sum([len(molecule_set) for molecule_set in self.molecules])

    def __getitem__(self, idx):
        return

