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


def test_buildgeoms_MoleculeSet():
    geom = pytest.mset[0]
    coords = geom.coords
    mol_labels = geom.labels
    atom_label_keys = geom.atoms[0].labels.keys()
    atom_labels = {key: tuple(atom.labels[key] for atom in geom.atoms) for key in atom_label_keys}
    new_geom = pytest.mset.build_geom(coords, mol_labels, atom_labels)
    assert geom.labels == new_geom.labels
    for atom1, atom2 in zip(geom.atoms, new_geom.atoms):
        assert atom1.at_num == atom2.at_num
        assert atom1.xyz == atom2.xyz
        assert atom1.labels == atom2.labels



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
