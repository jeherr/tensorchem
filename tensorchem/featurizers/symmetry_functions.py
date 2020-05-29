"""
Module for building features for model inputs from the molecular coordinates and atomic numbers
"""

import torch
from .util import dist_matrix_dense, cos_cutoff, gaussian_embed


def get_sym_funcs(params, at_nums, coords, elements):
    dist = dist_matrix_dense(coords)
    # padded_mol_mask = torch.ne(at_nums, 0)
    radial_embed = get_radial_embed(dist, params['r_nought'], params['eta'], params['rad_cut'])
    scatter = torch.eq(at_nums.unsqueeze(-1), elements)
    radial_channels = radial_embed.unsqueeze(-2) * scatter.unsqueeze(-1)
    return radial_channels


def get_radial_embed(dist, r_nought, eta, rad_cut):
    """
    Embeds a set of distances into a Gaussian function basis

    Args:
        dist: tensor of distances to be embedded in the Gaussian basis
        r_nought: a vector of distances to offset the Gaussian centers
        eta: parameter for controlling the width of the Gaussian functions
        rad_cut: cutoff distance beyond which all features should be zero

    Returns:
        radial_embed: coefficients of the distances embedded in the Gaussian basis
    """
    cutoff = cos_cutoff(dist, rad_cut)
    radial_embed = gaussian_embed(dist, r_nought, eta)
    return radial_embed * cutoff.unsqueeze(-1)
