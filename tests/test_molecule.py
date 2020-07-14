import tensorchem
import pytest

from tensorchem.dataset.molecule import MoleculeSet, Geometry, Atom

pytest.mset = MoleculeSet()
pytest.mset.load('tests/data/h2o.mset')


# MoleculeSet tests
def test_geom_MoleculeSet():
    for geom in pytest.mset.geometries:
        assert type(geom) == tensorchem.dataset.molecule.Geometry


def test_len_MoleculeSet():
    assert pytest.mset.__len__() == 2


def test_n_atoms_MoleculeSet():
    assert pytest.mset.n_atoms == 3


def test_at_nums_MoleculeSet():
    assert pytest.mset.at_nums == (8, 1, 1)


def test_at_symbs_MoleculeSet():
    assert pytest.mset.at_symbs == ('O', 'H', 'H')


def test_isomer_MoleculeSet():
    assert pytest.mset.is_isomer(pytest.mset)


def test_getitem_MoleculeSet():
    assert pytest.mset.__getitem__(1).n_atoms == 3


def test_save_MoleculeSet():
    pytest.mset.save()


# Geometry tests
def test_Geometry():
    for geom in pytest.mset.geometries:
        assert geom.n_atoms != 0 and type(geom.labels) == dict


def test_rep_Geometry():
    for geom in pytest.mset.geometries:
        assert type(geom.__repr__()) == str


def test_export_json_Geometry():
    for geom in pytest.mset.geometries:
        assert list(geom.export_json().keys()) == ['atoms', 'labels']


# Atom tests
def test_xyz_Atom():
    pytest.mset.atoms = [atom for geom in pytest.mset.geometries for atom in geom.atoms]
    for atom in pytest.mset.atoms:
        assert atom.xyz == (atom.x, atom.y, atom.z)


def test_export_json_Atom():
    pytest.mset.atoms = [atom for geom in pytest.mset.geometries for atom in geom.atoms]
    for atom in pytest.mset.atoms:
        assert list(atom.export_json().keys()) == ['atomic_num', 'xyz', 'labels']
