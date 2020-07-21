import pytest
import os
import json
from tensorchem.molecules import *

mol = Molecule.from_json(json.load(open(os.path.join(os.getcwd(),'tests/data/h2o.mset'), "r")))


# Molecule tests
def test_geom_Molecule():
    for geom in mol.geometries:
        assert type(geom) == Geometry


def test_len_Molecule():
    assert mol.__len__() == 1


def test_n_atoms_Molecule():
    assert mol.n_atoms == 3


def test_n_heavy_atoms_Molecule():
    assert mol.n_heavy_atoms == 1


def test_at_nums_Molecule():
    assert mol.atoms.tolist() == [8, 1, 1]


def test_at_symbs_Molecule():
    assert mol.at_symbs == ['O', 'H', 'H']


def test_getitem_Molecule():
    geom = mol[0]
    assert geom.n_atoms == mol.n_atoms
    assert geom.n_heavy_atoms == mol.n_heavy_atoms
    assert geom.at_symbs == mol.at_symbs
    assert set(geom.elements.tolist()) == set(mol.elements.tolist())


def test_save_Molecule():
    mol.save()


# Geometry tests
def test_Geometry():
    for geom in mol.geometries:
        assert geom.n_atoms != 0 and type(geom.labels) == dict


def test_rep_Geometry():
    for geom in mol.geometries:
        assert type(repr(geom)) == str


def test_export_json_Geometry():
    for geom in mol.geometries:
        assert set(geom.to_json().keys()) == {'atoms', 'xyz', 'labels'}
