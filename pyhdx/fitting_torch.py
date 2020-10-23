import torch.nn as nn
import torch.nn.functional as F
from torch.optim import SGD
import torch as t
from scipy import constants
import numpy as np
from pyhdx.models import Protein


class DeltaGModule(nn.Module):
    def __init__(self, deltaG):
        super(DeltaGModule, self).__init__()
        self.deltaG = deltaG

    def forward(self, inputs):
        temperature, X, k_int, timepoints = inputs

        pfact = t.exp(self.deltaG / (constants.R * temperature))
        uptake = 1 - t.exp(-t.matmul((k_int / (1 + pfact)), timepoints))
        return t.matmul(X, uptake)


class DeltaGFit(nn.Module):
    def __init__(self, deltaG):
        super(DeltaGFit, self).__init__()
        self.deltaG = deltaG

    def forward(self, inputs):
        """
        # inputs, list of:
            temperatures: scalar (1,)
            X (N_peptides, N_residues)
            k_int: (N_peptides, 1)

        """
        temperature, X, k_int, timepoints = inputs

        pfact = t.exp(self.deltaG / (constants.R * temperature))
        uptake = 1 - t.exp(-t.matmul((k_int / (1 + pfact)), timepoints))
        return t.matmul(X, uptake)

class TorchFitResult(object):
    def __init__(self, series, model, temperature=None, **metadata):
        self.series = series
        self.model = model
        self.temperature = temperature
        self.metadata = metadata

    @property
    def output(self):
        out_dict = {}
        out_dict['r_number'] = self.series.cov.r_number
        out_dict['deltaG'] = self.model.deltaG.detach().numpy().squeeze()
        if self.temperature is not None:
            pfact = np.exp(out_dict['deltaG'] / (constants.R * self.temperature))
            out_dict['pfact'] = pfact

        return Protein(out_dict, index='r_number')



