"""
Dataset objects for holding molecular geometries (atomic numbers, coordinates)
and their properties (energy, atomic forces, dipole moments) for training models.
"""

import json

import numpy as np
import torch

from torch.utils.data import Dataset as TorchDataset
from tensorchem.molecules import Molecule


class Dataset(TorchDataset):
    def __init__(self):
        super(Dataset, self).__init__()
        self.unique_atoms = []
        self.filename = None

    def _init_dataset(self):
        return

    # TODO develop a more careful save and load routine for JSON data
    def save(self, filename=None):
        # Will do some stuff to collect all necessary data into a format for JSON
        json_data = {
            "unique_atoms": self.unique_atoms
        }
        with open(filename, "w") as f:
            json.dump(json_data, f)

    def load(self, filename=None):
        with open(filename) as f:
            json_data = json.load(f)
        self.__dict__.update(json_data)


class MolDataset(Dataset):
    def __init__(self):
        super(MolDataset, self).__init__()
        self.samples = []  # Immutable type so order of molecules cannot change during training
        self.idx_map = {}  # Maps an overall sample index to the molecule and geometry indices

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
            self.samples.append(sample)

    @classmethod
    def from_mset(cls, msets):
        if type(msets) is Molecule:
            msets = [msets]
        mol_data = cls()
        for mset in msets:
            atoms = {"atomic_num": [atom.at_num for atom in mset.atoms],
                     "xyz": [atom.xyz for atom in mset.atoms]}
            atoms.update({key: value for atom in mset.atoms for key, value in atom.labels.items()})
            geometries = []
            for geom in mset.geometries:
                geoms = {"atomic_num": [atom.at_num for atom in geom.atoms],
                         "xyz": [atom.xyz for atom in geom.atoms]}
                geoms.update({key: value for key, value in geom.labels.items()})
                geometries.append(geoms)
            sample = {"atoms": atoms, "geometries": (geometries)} #geometries immutable
            mol_data.samples.append(sample)
        return mol_data


class MixedDataset(Dataset):
    def __init__(self):
        super(MixedDataset, self).__init__()
        self.samples = []
        self.filename = None

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
            self.samples.append(sample)

    @classmethod
    def from_mset(cls, msets):
        if type(msets) is Molecule:
            msets = [msets]
        mixed_data = cls()
        for mset in msets:
            for geom in mset.geometries:
                atoms = {"atomic_num": [atom.at_num for atom in geom.atoms],
                         "xyz": [atom.xyz for atom in geom.atoms]}
                atoms.update({key: value for atom in geom.atoms for key, value in atom.labels.items()})
                sample = {"atoms": atoms}
                sample.update({key: value for key, value in geom.labels.items()})
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
