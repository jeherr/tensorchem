"""
tensorchem is a Neural Network force field based on the High-Dimensional Neural Network Potential model developed
by Behler and Parrinello.
"""

import torch
import torch.nn as nn
import torch.optim as opt


class TensorChemTrainer(nn.Module):
    """
    Trainer for TensorChem networks. Coordinates the featurization of data to feed into TensorChem, recovers the latent
    featurization after the final non-linear layer to feed to linear regression layers for final label prediction.
    """

    def __init__(self, dataset, hyper_params):
        super().__init__()
        cuda_condition = torch.cuda.is_available()
        self.device = torch.device("cuda" if cuda_condition else "cpu")
        self.featurization = None  # TODO Need featurization modules implemented
        self.model = TensorChem(elements, layers, input_size)
        self.model.to(self.device)
        self.learn_rate = hyper_params['learning_rate']
        self.betas = hyper_params['betas']
        self.weight_decay = hyper_params['weight_decay']
        self.optimizer = opt.Adam(self.model.parameters(), lr=self.learn_rate, betas=self.betas,
                                  weight_decay=self.weight_decay)
        return


class TensorChem(nn.Module):
    """
    The network class which contains the necessary submodules

    Args:
        elements: a list of unique elements by atomic number in the training and testing data
        layers: list of ints defining the number of layers and hidden neurons per subnetwork
    """

    def __init__(self, elements, layers, input_size):
        super().__init__()
        subnets = [SubNet(element, input_size, layers) for element in elements]
        return

    def forward(self, at_num, coords):
        return


class SubNet(nn.Module):
    """
    Module for the element subnets.

    Args:
        layers: list of ints

    Returns:
        A feed-forward neural network. The number of hidden layers and neurons per layer is based on len(layers)
        and the values in the list, respectively.
    """

    def __init__(self, atomic_num, input_size, layers):
        super().__init__()
        self.layers = [nn.Linear(input_size, layers[0])]
        self.layers += [nn.Linear(layers[i], layers[i + 1]) for i in range(len(layers) - 1)]
        return

    def forward(self):
        return
