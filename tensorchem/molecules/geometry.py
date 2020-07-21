"""
Geometries are a single frame of a molecule with at least 3D coordinates (xyz)
and may also contain labels (e.g. potential energy, forces, dipoles)
"""
from typing import List

import torch

from qcelemental import periodictable
from qcelemental.models import Molecule as QCMol

FloatTensor = torch.FloatTensor
ByteTensor = torch.ByteTensor


class Geometry:
    def __init__(self, atoms: ByteTensor = None, xyz: FloatTensor = None, labels: dict = None):
        self.atoms = atoms
        self.xyz = xyz
        if labels is None:
            labels = {}
        self.labels = labels

    def __repr__(self) -> str:
        rep = f"{self.n_atoms}\n"
        for atom, xyz in zip(self.at_symbs, self.xyz.tolist()):
            rep += f"\n{atom}   {xyz[0]:9.6f}   {xyz[1]:9.6f}   {xyz[2]:9.6f}"
        return rep

    @property
    def at_symbs(self) -> List[str]:
        return [periodictable.to_symbol(atom) for atom in self.atoms.tolist()]

    @property
    def elements(self) -> set:
        return torch.unique(self.atoms)

    @property
    def n_atoms(self) -> int:
        return len(self.atoms.tolist())

    @property
    def n_heavy_atoms(self) -> int:
        return len([atom for atom in self.atoms.tolist() if atom != 1])

    @classmethod
    def from_json(cls, json_data: dict) -> 'Geometry':
        new_geom = cls()
        new_geom.atoms = torch.tensor(json_data["atoms"], dtype=torch.uint8)
        new_geom.xyz = torch.tensor(json_data["xyz"], dtype=torch.float32)
        new_geom.labels = {key: torch.tensor(label, dtype=torch.float32) for key, label in json_data["labels"].items()}
        return new_geom

    def to_json(self):
        data_dict = {"atoms": self.atoms.tolist(),
                     "xyz": self.xyz.tolist(),
                     "labels": {key: label.tolist() for key, label in self.labels.items()}
                     }
        return data_dict

    @classmethod
    def from_qcmol(cls, qcmol: QCMol) -> 'Geometry':
        new_geom = cls()
        new_geom.atoms = torch.tensor([periodictable.to_atomic_number(symbol) for symbol in qcmol.symbols])
        new_geom.xyz = torch.tensor(qcmol.geometry, dtype=torch.float32)
        return new_geom

    def write_xyz(self, filename: str) -> None:
        with open(filename, "w") as f:
            f.write(self.__repr__())
