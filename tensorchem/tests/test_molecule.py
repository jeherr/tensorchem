import tensorchem
import pytest
import sys
import json

from tensorchem.dataset.molecule import MoleculeSet, Geometry


def test_len_MoleculeSet():
    mset = MoleculeSet()
    mset.load('tensorchem/tests/h2o.mset')
    assert mset.__len__() == 2


def test_save_nofile_MoleculeSet():
    mset = MoleculeSet()
    mset.load('tensorchem/tests/h2o.mset')
    with pytest.raises(FileNotFoundError):
        mset.save()


def test_Geometry():
    with open('tensorchem/tests/h2o.mset', 'r') as data:
        data = json.loads(data.read())
    geom = Geometry(data['coordinates'], data['properties'])
    assert len(geom.coords) != 0 and type(geom.properties) == dict
