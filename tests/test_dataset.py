import tensorchem
import pytest

from tensorchem.dataset.dataset import MixedDataset
from tensorchem.molecules import Molecule


def test_len_MixedDataset():
    mixed_data = MixedDataset()
    mixed_data.load('tests/data/h2o.dset')
    assert mixed_data.__len__() == 2


def test_getitem_MixedDataset():
    mixed_data = MixedDataset()
    mixed_data.load('tests/data/h2o.dset')
    assert list(mixed_data.__getitem__(0).keys()) == ["atomic_numbers", "coordinates", "wb97x-d.6-311gss.mulliken_charge"]


def test_save_nofile_MixedDataset():
    mixed_data = MixedDataset()
    mixed_data.load('tests/data/h2o.dset')
    with pytest.raises(FileNotFoundError):
        mixed_data.save()


#def test_MixedDataset_from_mset():
#    mset = Molecule()
#    mset.load('h2o.mset', './tests/data')
#    assert type(MixedDataset.from_mset(mset)) == tensorchem.dataset.dataset.MixedDataset
