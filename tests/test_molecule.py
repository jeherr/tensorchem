import tensorchem
import pytest
import sys
import json

from tensorchem.dataset.molecule import MoleculeSet, Geometry, Atom

pytest.mset = MoleculeSet()
pytest.mset.load('tests/data/h2o.mset')

def test_geom_MoleculeSet():
    for geom in pytest.mset.geometries:
        assert type(geom) == tensorchem.dataset.molecule.Geometry

def test_len_MoleculeSet():
    assert pytest.mset.__len__() == 2

def test_getitem_MoleculeSet():
    assert pytest.mset.__getitem__(1).n_atoms == 3

def test_eq_MoleculeSet():
    mset1, mset2 = MoleculeSet(), MoleculeSet()
    mset1.load('tests/data/h2o.mset')
    mset2.load('tests/data/h2o.mset')
    assert mset1 == mset2

def test_save_nofile_MoleculeSet():
    with pytest.raises(FileNotFoundError):
        pytest.mset.save()

def test_Geometry():
    for geom in pytest.mset.geometries:
        assert geom.n_atoms != 0 and type(geom.labels) == dict

def test_rep_Geometry():
    for geom in pytest.mset.geometries:
        assert type(geom.__repr__()) == str

def test_export_json_Geometry():
    for geom in pytest.mset.geometries:
        assert list(geom.export_json().keys()) == ['atoms', 'labels']

def test_xyz_Atom():
    pytest.mset.atoms = [atom for geom in pytest.mset.geometries for atom in geom.atoms]
    for atom in pytest.mset.atoms:
        assert atom.xyz == [atom.x, atom.y, atom.z]

def test_export_json_Atom():
    pytest.mset.atoms = [atom for geom in pytest.mset.geometries for atom in geom.atoms]
    for atom in pytest.mset.atoms:
        assert list(atom.export_json().keys()) == ['atomic_num', 'xyz', 'labels'] 
