"""
Functions for generating samples of molecules, such as Normal Modes Sampling
"""

import torch
from qcelemental import periodictable, physical_constants


def sample_normal_modes(conformer, n_samples, temperature):
    """
    Take a conformer and return randomly sampled geometries by applying
    random displacements according to the normal modes of the molecule
    """
    atoms = conformer.atoms
    n_atoms = atoms.shape[0]
    original_xyz = conformer.xyz
    masses = [periodictable.to_mass(atom) for atom in atoms.tolist()]
    normal_modes = conformer.labels['normal modes']

    force_constants = conformer.labels['force constants']
    n_modes = len(normal_modes)
    kb = physical_constants.pc['boltzmann constant']

    new_coords = []
    for sample in range(n_samples):
        cis = torch.rand(n_modes)
        c_scale = torch.rand(1)
        cis = c_scale * cis / cis.sum()
        displacement_scalars = (((torch.bernoulli(torch.fill(n_modes, 0.5)) - 0.5) * 2.0)
                               * torch.sqrt(3.0 * cis * n_atoms * kb * temperature / force_constants))
        displacements = [mode * displacement_scalars[i] for i, mode in enumerate(normal_modes)]


def apply_transformation(coords, transformation):
    transformation.mean(dim=1)
