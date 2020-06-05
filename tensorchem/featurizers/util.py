"""
Functions common to several featurizers
"""

import torch
from numpy import pi


def cos_cutoff(dist, cutoff):
    """
    Cosine cutoff function used to smoothly decay feature signals to 0 at the limit of the sensory range of the
    features

    Args:
          dist: Array of distances between two points in angstroms
          cutoff: distance beyond which all values should go to 0

    Returns:
          scaling_factor: Array with the same shape as dist. Used to scale features based on their distance
    """
    cos_factor = 0.5 * (torch.cos(pi * dist / cutoff) + 1.0)
    cos_factor = torch.where(dist < cutoff, cos_factor, torch.zeros_like(cos_factor))
    return cos_factor


def cos_angle(dxyz_ij, dxyz_ik):
    """
    Returns the angle between two vectors with the same initial point

    Args:
        dxyz_ij: Array of three dimensional vectors between the initial points and the first terminal point
        dxyz_ik: Array of three dimensional vectors between the initial points and the second terminal point

    Returns:
         angle: Array of angles between the ij vector and ik vector
    """
    dist_ij, dist_ik = torch.norm(dxyz_ij, dim=-1), torch.norm(dxyz_ik, dim=-1)
    dist_ij_ik = dist_ij * dist_ik
    ij_dot_ik = torch.sum(dxyz_ij * dxyz_ik, dim=-1)
    cos_ij_ik = ij_dot_ik / dist_ij_ik
    return torch.acos(cos_ij_ik)


def gaussian_embed(dist, gauss_offset, gauss_width):
    """
    Embeds the distance between two points into a gaussian vector basis.

    Args:
        dist: Array of distances between two points in angstroms
        gauss_offset: Array of distance offsets for the center of the Gaussian embedding functions
        gauss_width: Array of Gaussian width parameters for the embedding functions. Must have the same shape as
            gauss_offset

    Returns:
        gauss_embeds: Array with the same leading dimensions as dist but has an extra dimension for the number of
            Gaussian peaks the distance is being embedded in.
    """
    exponent = -1.0 * gauss_width * torch.square(dist.unsqueeze(-1) - gauss_offset)
    return torch.exp(exponent)


def dist_matrix_dense(coords):
    """
    Calculates the distance matrix for a batch of atomic coordinates matrices

    Args:
        coords: ... x Na x 3 torch tensor of atomic positions

    Returns:
        dist: ... x Na x Na torch tensor of distances between atoms
    """
    d_xyz = coords.unsqueeze(-2) - coords.unsqueeze(-3)
    return torch.norm(d_xyz, dim=-1)
