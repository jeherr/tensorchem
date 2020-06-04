"""
Dataset objects for holding molecular geometries (atomic numbers, coordinates)
and their labels (energy, atomic forces, dipole moments).
"""

import json
from typing import List, Tuple, Optional

from ase.data import chemical_symbols


class MoleculeSet:
    def __init__(self, atoms: Tuple['Atom', ...] = None):
        self.atoms = atoms
        self.trajectories = {}
        self.identifiers = {}
        self.min_geom = None
        self.filename = None

    def __len__(self) -> int:
        return len(self.geometries)

    def __getitem__(self, idx: int) -> List[List[int]]:
        return self.geometries[idx]

    def __hash__(self) -> hash(tuple):
        return hash(tuple([atom.at_num for atom in self.atoms]))

    def compare_hash(self, other: 'MoleculeSet') -> bool:
        # TODO We should probably rename this to a more descriptive named method instead of __eq__. It might be
        #  confusing to users. Equality comparison is different from what we are doing here, so it may not be very
        #  "Pythonic"
        if isinstance(other, MoleculeSet):
            return self.__hash__() == other.__hash__()
        return NotImplemented

    @property
    def geometries(self) -> Tuple[List[List[int]]]:
        return list(self.trajectories.values())[0]

    @property
    def at_nums(self) -> Tuple[int, ...]:
        return tuple([atom.at_num for atom in self.atoms])

    @property
    def at_symbs(self) -> Tuple[str, ...]:
        return tuple([chemical_symbols[atom.at_num] for atom in self.atoms])

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def min_geom(self) -> 'Geometry':
        energies = [value if "energy" in key else 0 for geom in self.geometries for key, value in geom.labels.items()]
        return self.geometries[index(min(energies))]

    def build_geom(self, coords: list, mol_labels: dict, atom_labels: dict) -> 'Geometry':
        geom_atoms = tuple([Atom(atom.at_num) for atom in self.atoms])
        for i, atom in enumerate(geom_atoms):
            atom.xyz = (coords[i][0], coords[i][1], coords[i][2])
        for key, label in atom_labels.items():
            for i, atom in enumerate(geom_atoms):
                atom.labels[key] = label[i]
        geom = Geometry(geom_atoms)
        geom.labels = mol_labels
        return geom

    def save(self, filename: str = None):
        if filename is None:
            if self.filename is None:
                raise FileNotFoundError("No filename given for saving")
            else:
                filename = self.filename
        json_data = {"atoms": [atom.export_json() for atom in self.atoms],
                     "trajectories": {key: [geom.export_json() for geom in value] for key, value in self.trajectories.items()},
                     "identifiers": self.identifiers
        }
        with open(filename, 'w') as f:
            json.dump(json_data, f)

    def load(self, filename: str = None):
        if filename is None:
            if self.filename is None:
                raise FileNotFoundError("No filename given for loading")
            else:
                filename = self.filename
        with open(filename) as f:
            json_data = json.load(f)
        self.atoms = tuple([Atom.from_json(data) for data in json_data['atoms']])
        self.identifiers = json_data['identifiers']
        for key, value in json_data['trajectories'].items():
            self.trajectories.update({key: tuple([Geometry.from_json(geom_data) for geom_data in value])})


class Geometry:
    def __init__(self, atoms: List['Atom'] = None, labels: dict = None):
        self.atoms = atoms
        if labels is None:
            labels = {}
        self.labels = labels

    def __repr__(self) -> str:
        # TODO add some of that fancy printing crap to make this print the xyz in a nicer format
        rep = str(self.n_atoms)
        rep += "\n\n"
        for atom in self.atoms:
            rep += "     ".join((str(atom.at_num), str(atom.x), str(atom.y), str(atom.z)))
            rep += "\n"
        return rep

    @property
    def at_nums(self) -> Tuple[Optional[int], ...]:
        return tuple([atom.at_num for atom in self.atoms])

    @property
    def at_symbs(self) -> Tuple[Optional[str], ...]:
        return tuple([chemical_symbols[atom.at_num] for atom in self.atoms])

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def coords(self) -> Tuple[Optional[Tuple[float, float, float]], ...]:
        return tuple([atom.xyz for atom in self.atoms])

    @classmethod
    def from_json(cls, json_data: dict) -> 'Geometry':
        new_geom = cls()
        new_geom.atoms = tuple([Atom.from_json(atom_data) for atom_data in json_data['atoms']])
        new_geom.labels = json_data['labels']
        return new_geom

    def export_json(self) -> dict:
        return {"atoms": [atom.export_json() for atom in self.atoms],
                "labels": self.labels}


class Atom:
    def __init__(self, at_num: int = None, xyz: Tuple[float, float, float] = (None, None, None)):
        self.at_num = at_num
        self.xyz = xyz
        self.labels = {}

    @property
    def x(self) -> float:
        return self.xyz[0]

    @property
    def y(self) -> float:
        return self.xyz[1]

    @property
    def z(self) -> float:
        return self.xyz[2]

    @classmethod
    def from_json(cls, json_data: dict) -> 'Atom':
        new_atom = cls()
        new_atom.at_num = json_data['atomic_num']
        new_atom.xyz = tuple([json_data['xyz']])
        new_atom.labels = json_data['labels']
        return new_atom

    def export_json(self) -> dict:
        return {"atomic_num": self.at_num,
                "xyz": self.xyz,
                "labels": self.labels}
