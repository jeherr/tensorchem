"""
Dataset objects for holding molecular geometries (atomic numbers, coordinates)
and their properties (energy, atomic forces, dipole moments) for training models.
"""

import json

import numpy as np
import torch

from torch.utils.data import Dataset as TorchDataset
from tensorchem.dataset.molecule import MoleculeSet, Geometry


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
        self.samples = ()  # Immutable type so order of molecules cannot change during training
        self.idx_map = {}  # Maps an overall sample index to the molecule and geometry indices

    def __len__(self):
        return sum([len(molecule_set) for molecule_set in self.samples])

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
            "properties": [geom.properties for sample in self.samples for geom in sample['geometries']]
        }
        with open(filename, "w") as f:
            json.dump(json_data, f)

    def load(self, filename=None):
        if filename is None:
            if self.filename is None:
                raise FileNotFoundError("No filename given for loading")
            else:
                filename = self.filename
        with open(filename, 'r') as f:
            for line in f:
                mset = json.loads(line)
                sample = {
                    "atomic_number": mset['atomic_number'],
                    "geometries": [Geometry(np.array(mset['coordinates']), mset['properties'])]
                }
                self.samples.append(sample)


class MixedDataset(Dataset):
    def __init__(self):
        super(MixedDataset, self).__init__()
        self.samples = []

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, idx):
        item = {}
        for key, value in self.samples[idx].items():
            if type(value) == np.ndarray:
                item.update({key: torch.from_numpy(value)})
            elif type(value) == float:
                item.update({key: torch.FloatTensor([value])})
            else:
                item.update({key: torch.FloatTensor(value)})
        return item

    def save(self, filename=None):
        if filename is None:
            if self.filename is None:
                raise FileNotFoundError("No filename given for saving")
            else:
                filename = self.filename
        json_data = []
        for sample in self.samples:
            data = {}
            for key, value in sample.items():
                if type(value) is np.ndarray:
                    data.update({key: value.tolist()})
                else:
                    data.update({key: value})
            json_data.append(data)
        with open(filename, "w") as f:
            json.dump(json_data, f)

    def load(self, filename=None):
        if filename is None:
            if self.filename is None:
                raise FileNotFoundError("No filename given for loading")
            else:
                filename = self.filename
        with open(filename, "r") as f:
            json_data = json.loads(f.read())
        if type(json_data) is not list:
            json_data = [json_data]
        for data in json_data:
            sample = {}
            for key, value in data.items():
                if type(value) == list:
                    sample.update({key: value})
                elif type(value) == dict:
                    for k, v in value.items():
                        sample.update({k: v})
            for key, value in sample.items():
                if key == "coordinates":
                    for data in value:
                        sample.update({key: np.array(data)})
                        self.samples.append(sample)

    @classmethod
    def from_mset(cls, msets):
        if type(msets) is MoleculeSet:
            msets = [msets]
        mixed_data = cls()
        for mset in msets:
            for geom in mset.geometries:
                sample = {"atomic_numbers": mset.atomic_nums, "coordinates": geom.coords}
                sample.update({key: value for key, value in geom.properties.items()})
                mixed_data.samples.append(sample)
        return mixed_data


class Sample:
    def __init__(self, atomic_nums, coords, labels):
        self.atomic_nums = torch.tensor(atomic_nums, dtype=torch.uint8)
        self.coords = torch.tensor(coords, dtype=torch.float32)
        self.labels = {}
        for key, value in labels.items():
            if type(value) == np.ndarray:
                self.labels.update({key: torch.from_numpy(value)})
            elif type(value) == float:
                self.labels.update({key: torch.FloatTensor([value])})
            else:
                self.labels.update({key: torch.FloatTensor(value)})



if __name__ == "__main__":
    mol1 = MoleculeSet()
    mol1.load('../data/ani1x-mol.mset')
    ds1 = MixedDataset.from_mset(mol1)
    print(ds1.__len__())
