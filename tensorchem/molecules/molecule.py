import os
import json
from typing import List

import torch

from qcelemental import periodictable

from tensorchem.molecules.geometry import Geometry
from tensorchem.molecules.trajectory import Trajectory

ByteTensor = torch.ByteTensor


class Molecule:
    def __init__(self, atoms: ByteTensor = None, charge: int = None, identifiers: dict = None,
                 trajectories: list = None, filename: str = None, filepath: str = None):
        self.atoms = atoms
        self.charge = charge
        if identifiers is None:
            identifiers = {}
        self.identifiers = identifiers
        if trajectories is None:
            trajectories = []
        self.trajectories = trajectories
        self.filename = filename
        self.filepath = filepath

    def __len__(self) -> int:
        return len(self.geometries)

    def __getitem__(self, idx: int) -> 'Geometry':
        return self.geometries[idx]

    @property
    def chem_formula(self):
        formula = {}
        for atom in self.atoms.tolist():
            if periodictable.to_symbol(atom) in formula.keys():
                formula[periodictable.to_symbol(atom.at_num)] += 1
            else:
                formula[periodictable.to_symbol(atom.at_num)] = 1
        return formula

    @property
    def at_symbs(self) -> List[str]:
        return [periodictable.to_symbol(atom) for atom in self.atoms.tolist()]

    @property
    def elements(self) -> set:
        return torch.unique(self.atoms)

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def n_heavy_atoms(self) -> int:
        return len([atom for atom in self.atoms.tolist() if atom != 1])

    @property
    def geometries(self) -> List['Geometry']:
        return [geom for trajectory in self.trajectories for geom in trajectory.geometries]

    def to_json(self):
        json_data = {"atoms": self.atoms.tolist(),
                     "charge": self.charge,
                     "trajectories": [trajectory.to_json() for trajectory in self.trajectories],
                     "identifiers": self.identifiers,
                     "filename": self.filename,
                     "filepath": self.filepath
                     }
        return json_data

    def save(self, filename: str = None, filepath: str = None):
        if filename is None:
            if self.filename is None:
                raise FileNotFoundError("No filename given for saving")
            else:
                filename = self.filename
        if filepath is None:
            if self.filepath is None:
                raise FileNotFoundError("No filepath given for saving")
            else:
                filepath = self.filepath
        json_data = self.to_json()
        with open(os.path.join(filepath, filename), 'w') as f:
            json.dump(json_data, f)

    @classmethod
    def from_json(cls, json_data):
        new_mol = cls()
        new_mol.atoms = torch.tensor(json_data["atoms"], dtype=torch.uint8)
        new_mol.charge = json_data["charge"]
        new_mol.identifiers = json_data["identifiers"]
        new_mol.trajectories = [Trajectory.from_json(traj_data) for traj_data in json_data["trajectories"]]
        new_mol.filename = json_data["filename"]
        new_mol.filepath = json_data["filepath"]
        return new_mol

    def load(self, filename: str = None, filepath: str = None):
        if filename is None:
            if self.filename is None:
                raise FileNotFoundError("No filename given for loading")
            else:
                filename = self.filename
        if filepath is None:
            if self.filepath is None:
                raise FileNotFoundError("No filepath given for loading")
            else:
                filepath = self.filepath
        with open(os.path.join(filepath, filename), "r") as f:
            json_data = json.load(f)
            new_mol = self.from_json(json_data)
            self.__dict__.update(new_mol.__dict__)
        return

    @classmethod
    def from_qcschema(cls, qcschema):
        new_mol = cls()
        topology = qcschema["molecule"]
        at_symbs = topology["symbols"]
        new_mol.atoms = [periodictable.to_atomic_number(atom) for atom in at_symbs]
        new_mol.identifiers = {}
        new_mol.trajectories = []
        new_mol.charge = topology["molecular_charge"] if "molecular_charge" in topology.keys() else None
        new_mol.bonds = topology["connectivity"] if "connectivity" in topology.keys() else []
        new_mol.multiplicity = topology["molecular_multiplicity"] if "molecular_multiplicity" in topology.keys() \
            else None
        if "name" in topology.keys():
            new_mol.identifiers["name"] = [topology["name"]]
        return new_mol
