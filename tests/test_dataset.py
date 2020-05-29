import tensorchem
import pytest
import sys
import numpy as np

from tensorchem.dataset.dataset import MixedDataset
from tensorchem.dataset.molecule import MoleculeSet


def test_len_MixedDataset():
    mixed_data = MixedDataset()
    mixed_data.load('tests/h2o.dset')
    assert mixed_data.__len__() == 2


def test_getitem_MixedDataset():
    mixed_data = MixedDataset()
    mixed_data.load('tests/h2o.dset')
    assert list(mixed_data.__getitem__(0).keys()) == ["atomic_numbers", "coordinates", "deg_of_freedom"]


def test_save_nofile_MixedDataset():
    mixed_data = MixedDataset()
    mixed_data.load('tests/h2o.dset')
    with pytest.raises(FileNotFoundError):
        mixed_data.save()


def test_MixedDataset_from_mset():
    mset = MoleculeSet()
    mset.load('tests/h2o.mset')
    assert type(MixedDataset.from_mset(mset)) == tensorchem.dataset.dataset.MixedDataset