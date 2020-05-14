"""
Thanks to the following citation and links within to available code for this downloaded
Smith, J.S., Zubatyuk, R., Nebgen, B. et al. The ANI-1ccx and ANI-1x data sets, coupled-cluster and density
functional theory properties for molecules. Sci Data 7, 134 (2020). https://doi.org/10.1038/s41597-020-0473-z
"""

import h5py
import numpy as np
import json
import sys

# Path to the ANI-1x data set
path_to_h5file = '/home/adriscoll/tensormol-jax/tensormol_jax/data/ani1x-release.h5'

argumentList = sys.argv
argumentList.pop(0)

# List of keys to point to requested data
data_keys = np.array(argumentList)

def iter_data_buckets(h5filename, keys=data_keys):
    """ Iterate over buckets of data in ANI HDF5 file.
    Yields dicts with atomic numbers (shape [Na,]) coordinated (shape [Nc, Na, 3])
    and other available properties specified by `keys` list, w/o NaN values.
    """
    keys_set = set(keys)
    keys_set.discard('atomic_numbers')
    keys_set.discard('coordinates')
    with h5py.File(h5filename, 'r') as f:
        for grp in f.values():
            Nc = grp['coordinates'].shape[0]
            mask = np.ones(Nc, dtype=np.bool)
            data = dict((k, grp[k][()]) for k in keys)
            for k in keys:
                v = data[k].reshape(Nc, -1)
                mask = mask & ~np.isnan(v).any(axis=1)
            if not np.sum(mask):
                continue
            d = dict((k, data[k][mask]) for k in keys)
            d['atomic_numbers'] = grp['atomic_numbers'][()]
            d['coordinates'] = grp['coordinates'][()][mask]
            yield d


json_data = []
molecule_dict = {}

name_dict = {
    'coordinates':"atomic_positions",
    'atomic_numbers':"atomic_numbers",
    'wb97x_dz.energy':"total_dz_energy",
    'wb97x_tz.energy':"total_tz_energy",
    'ccsd(t)_cbs.energy':"total_cbs_energy",
    'hf_dz.energy':"hf_dz_energy",
    'hf_tz.energy':"hf_tz_energy",
    'hf_qz.energy':"hf_qz_energy",
    'npno_ccsd(t)_dz.corr_energy':"npno_dz_energy",
    'npno_ccsd(t)_tz.corr_energy':"npno_tz_energy",
    'npno_ccsd(t)_qz.corr_energy':"npno_qz_energy",
    'mp2_dz.corr_energy':"mp2_dz_energy",
    'mp2_tz.corr_energy':"mp2_tz_energy",
    'mp2_qz.corr_energy':"mp2_qz_energy",
    'wb97x_dz.forces':"atomic_dz_forces",
    'wb97x_tz.forces':"atomic_tz_forces",
    'wb97x_dz.dipole':"dipole_dz",
    'wb97x_tz.dipole':"dipole_tz",
    'wb97x_tz.quadrupole':"moments",
    'wb97x_dz.cm5_charges':"cm5_charges",
    'wb97x_dz.hirshfeld_charges':"hirshfeld_charges",
    'wb97x_dz.mbis_charges':"mbis_charges",
    'wb97x_dz.mbis_dipoles':"mbis_dipoles",
    'wb97x_dz.mbis_quadrupoles':"mbis_electric",
    'wb97x_dz.mbis_octupoles':"mbis_moments",
    'wb97x_dz.mbis_volumes':"atomic_volumes"
}
    

# Example for extracting DFT/DZ energies and forces
for data in iter_data_buckets(path_to_h5file, keys=data_keys):
    X = data['coordinates']
    Y = data['atomic_numbers']
    molecule_dict.update({
        "atomic_positions":X.tolist(),
        "atomic_numbers":Y.tolist()
    }) 
    for key in data_keys:
        if (key == 'coordinates') or (key == 'atomic_positions'):
            continue
        molecule_dict.update({
            name_dict[key] : data[key].tolist()
        }) 
    json_data.append(molecule_dict)

with open('/home/adriscoll/tensormol-jax/tensormol_jax/data/ani1x-release.json', 'w') as outfile:
    json.dump(json_data, outfile, ensure_ascii=False)

