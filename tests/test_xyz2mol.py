import tensorchem
import pytest

from tensorchem.util.xyz2mol import geom_to_smiles
from tensorchem.dataset.molecule import MoleculeSet

pytest.mset = MoleculeSet()
pytest.mset.load('tests/data/h2o.mset')


def test_geom_to_smiles():
    geom = pytest.mset.geometries[0]
    smiles = geom_to_smiles(geom)
    assert smiles == "O"
