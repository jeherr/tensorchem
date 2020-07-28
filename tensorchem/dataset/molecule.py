"""
Dataset objects for holding molecular geometries (atomic numbers, coordinates)
and their labels (energy, atomic forces, dipole moments).
"""

import json
from typing import List, Tuple, Optional

from ase.data import chemical_symbols

from .labels import Potential, Forces, Dipole, Quadrupole, Charge


class MoleculeSet:
    def __init__(self, atoms: Tuple['Atom', ...] = None):
        self.atoms = atoms
        self.formal_charge = None
        self.trajectories = {}
        self.identifiers = {}
        self.min_geom = None
        self.filename = None

    def __len__(self) -> int:
        return len(self.geometries)

    def __getitem__(self, idx: int) -> 'Geometry':
        return self.geometries[idx]

    def is_isomer(self, other: 'MoleculeSet') -> bool:
        if isinstance(other, MoleculeSet):
            return self.chem_formula == other.chem_formula
        return NotImplemented

    @property
    def chem_formula(self):
        formula = {}
        for atom in self.atoms:
            if chemical_symbols[atom.at_num] in formula.keys():
                formula[chemical_symbols[atom.at_num]] += 1
            else:
                formula[chemical_symbols[atom.at_num]] = 1
        return formula

    @property
    def geometries(self) -> List['Geometry']:
        return [geom for value in self.trajectories.values() for geom in value]

    @property
    def at_nums(self) -> Tuple[int, ...]:
        return tuple([atom.at_num for atom in self.atoms])

    @property
    def at_symbs(self) -> Tuple[str, ...]:
        return tuple([chemical_symbols[atom.at_num] for atom in self.atoms])

    @property
    def elements(self) -> set:
        return set([self.at_symbs])

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

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
        if self.filename is None:
            self.filename = filename
        with open(filename) as f:
            json_data = json.load(f)
        self.atoms = tuple([Atom.from_json(data) for data in json_data['atoms']])
        self.identifiers = json_data['identifiers']
        for key, value in json_data['trajectories'].items():
            self.trajectories.update({key: tuple([Geometry.from_json(geom_data) for geom_data in value])})

    def write_xyz_trajectory(self):
        with open(".".join((self.filename.rstrip(".mset"), "xyz")), "w") as f:
            for geom in self.geometries:
                f.write(str(geom))
                f.write("\n")


class Geometry:
    def __init__(self, atoms: List['Atom'] = None, labels: dict = None):
        self.atoms = atoms
        self.labels = {}
        if labels is not None:
            self.add_labels(labels)

    def __repr__(self) -> str:
        # TODO add some of that fancy printing crap to make this print the xyz in a nicer format
        rep = f"{self.n_atoms}\n"
        for atom in self.atoms:
            rep += f"\n{atom.at_symb}   {atom.x:9.6f}   {atom.y:9.6f}   {atom.z:9.6f}"
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
        new_geom.add_labels(json_data['labels'])
        return new_geom

    def add_labels(self, mol_labels: list):
        for key, value in mol_labels.items():
            description = key.split('.')
            label_type = description[0]
            functional = description[-2]
            basis = description[-1]
            if label_type == "potential":
                value = float(value)
                mol_label = Potential(value, functional, basis)
            elif label_type == "dipole":
                value = tuple(value)
                mol_label = Dipole(value, functional, basis)
            elif label_type == "quadrupole":
                value = tuple(value)
                mol_label = Quadrupole(value, functional, basis)
            else:
                raise RuntimeError(f"No label type know for {label_type}")
            if label_type in self.labels.keys():
                self.labels[label_type].append(mol_label)
            else:
                self.labels[label_type] = [mol_label]

    def export_json(self) -> dict:
        labels_dict = {}
        for key, value in self.labels.items():
            for label in value:
                ke, val = label.export_json()
                labels_dict[ke] = val
        return {"atoms": [atom.export_json() for atom in self.atoms],
                "labels": labels_dict}

    def write_xyz(self, filename: str) -> None:
        with open(filename, "w") as f:
            f.write(self.__repr__())


class Atom:
    def __init__(self, at_num: int = None, xyz: Tuple[float, float, float] = (None, None, None)):
        self.at_num = at_num
        self.xyz = xyz
        self.labels = {}

    @property
    def at_symb(self) -> str:
        return chemical_symbols[self.at_num]

    def __repr__(self) -> str:
        return self.at_symb

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
        if json_data['xyz'] is None:
            new_atom.xyz = (None, None, None)
        else:
            new_atom.xyz = tuple(json_data['xyz'])
        new_atom.add_labels(json_data['labels'])
        return new_atom

    def add_labels(self, atom_labels: dict):
        for key, value in atom_labels.items():
            if key == 'atomic_num' or key == 'xyz':
                continue
            description = key.split('.')
            label_type = description[0]
            functional = description[-2]
            basis = description[-1]
            if label_type == "forces":
                value = tuple(value)
                atom_label = Forces(value, functional, basis)
            elif label_type == "charge":
                value = float(value)
                partitioning = description[1]
                atom_label = Charge(value, partitioning, functional, basis)
            else:
                raise RuntimeError(f"No label type know for {label_type}")
            if label_type in self.labels.keys():
                self.labels[label_type].append(atom_label)
            else:
                self.labels[label_type] = [atom_label]

    def export_json(self) -> dict:
        labels_dict = {}
        for key, value in self.labels.items():
            for label in value:
                ke, val = label.export_json()
                labels_dict[ke] = val
        return {"atomic_num": self.at_num,
                "xyz": self.xyz,
                "labels": labels_dict}
