import pytest

from tensorchem.featurizers.util import *


def test_cos_cutoff():
    scaling = cos_cutoff(torch.arange(0.0, 7.5, step=0.5), torch.tensor(5.0, dtype=torch.float32))
    assert torch.allclose(scaling, torch.tensor([1.0000, 0.975528, 0.904508, 0.793893, 0.654508, 0.5000, 0.345492,
                                                 0.206107, 0.0954915, 0.0244717, 0.0, 0.0, 0.0, 0.0, 0.0]))


def test_cos_angle():
    dxyz1 = torch.tensor([[1, 1, 1], [1, 1, 1], [1, 1, 0], [1, 1, 1], [2, 2, 2], [-1, 1, 1]], dtype=torch.float32)
    dxyz2 = torch.tensor([[1, 1, 1], [1, 1, 0], [1, 1, 1], [1, 0, 0], [2, 2, 2], [1, 0, 0]], dtype=torch.float32)
    angles = cos_angle(dxyz1, dxyz2)
    assert torch.allclose(angles, torch.tensor([0.00000, 0.61548, 0.61548, 0.95532, 0.00000, 2.18628],
                                               dtype=torch.float32))

