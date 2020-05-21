import tensormol_jax
import pytest
import sys
import numpy as np

from tensormol_jax.dataset.dataset import MixedDataset
from tensormol_jax.dataset.molecule import MoleculeSet

def test_len_MixedDataset():
    mixed_data = MixedDataset()
    mixed_data.load('tensormol_jax/tests/h2o.mset')
    assert mixed_data.__len__() == 1

def test_getitem_MixedDataset():
    mixed_data = MixedDataset()
    mixed_data.load('tensormol_jax/tests/h2o.mset')
    assert list(mixed_data.__getitem__(0).keys()) == ["atomic_number", "coordinates", "deg_of_freedom"]

def test_save_nofile_MixedDataset():
    mixed_data = MixedDataset()
    mixed_data.load('tensormol_jax/tests/h2o.mset')
    with pytest.raises(FileNotFoundError):
        mixed_data.save()

def test_MixedDataset_from_mset():
    MSet = MoleculeSet()
    MSet.load('tensormol_jax/tests/h2o.mset')
    assert type(MixedDataset.from_mset(MSet)) == tensormol_jax.dataset.dataset.MixedDataset
