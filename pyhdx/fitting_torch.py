import torch.nn as nn
import torch.nn.functional as F
from torch.optim import SGD
import torch as t
from scipy import constants
import numpy as np
from pyhdx.models import Protein


class DeltaGFit(nn.Module):
    def __init__(self, deltaG):
        super(DeltaGFit, self).__init__()
        self.deltaG = deltaG

    def forward(self, temperature, X, k_int, timepoints):
        """
        # inputs, list of:
            temperatures: scalar (1,)
            X (N_peptides, N_residues)
            k_int: (N_peptides, 1)

        """

        pfact = t.exp(self.deltaG / (constants.R * temperature))
        uptake = 1 - t.exp(-t.matmul((k_int / (1 + pfact)), timepoints))
        return t.matmul(X, uptake)


class TorchFitResult(object):
    def __init__(self, series, model, temperature=None, **metadata):
        self.series = series
        self.model = model
        self.temperature = temperature
        self.metadata = metadata

    def estimate_errors(self):
        # boolean array to select residues which are exchanging (ie no nterminal resiudes, no prolines, no regions without coverage)
        bools = self.series.cov['exchanges'].to_numpy()
        r_number = self.series.cov.r_number[bools]  # Residue number which exchange

        dtype = t.float64
        temperature = t.tensor([self.temperature], dtype=dtype)
        X = t.tensor(self.series.cov.X[:, bools], dtype=dtype)  # Np x Nr, non-exchanging residues removed
        k_int = t.tensor(self.series.cov['k_int'][bools].to_numpy(), dtype=dtype).unsqueeze(-1)  # Nr x 1
        timepoints = t.tensor(self.series.timepoints, dtype=dtype).unsqueeze(0)  # 1 x Nt

        deltaG = t.tensor(self.deltaG[bools], dtype=dtype)
        output_data = t.tensor(self.series.uptake_corrected.T, dtype=dtype)

        def calc_loss(deltaG_input):
            criterion = t.nn.MSELoss(reduction='sum')
            pfact = t.exp(deltaG_input.unsqueeze(-1) / (constants.R * temperature))
            uptake = 1 - t.exp(-t.matmul((k_int / (1 + pfact)), timepoints))
            output = t.matmul(X, uptake)

            loss = criterion(output, output_data)
            return loss

        hessian = t.autograd.functional.hessian(calc_loss, deltaG)
        hessian_inverse = t.inverse(-hessian)
        covariance = np.sqrt(np.abs(np.diagonal(hessian_inverse)))

        return Protein({'covariance': covariance, 'r_number': r_number}, index='r_number')

    @property
    def deltaG(self):
        return self.model.deltaG.detach().numpy().squeeze()

    @property
    def output(self):
        out_dict = {}
        out_dict['r_number'] = self.series.cov.r_number
        out_dict['sequence'] = self.series.cov['sequence'].to_numpy()
        out_dict['_deltaG'] = self.deltaG
        out_dict['deltaG'] = out_dict['_deltaG'].copy()
        out_dict['deltaG'][~self.series.cov['exchanges']] = np.nan
        if self.temperature is not None:
            pfact = np.exp(out_dict['deltaG'] / (constants.R * self.temperature))
            out_dict['pfact'] = pfact

        #todo add possibility to add append series to protein?
        #todo update order of columns
        protein = Protein(out_dict, index='r_number')
        protein_cov = self.estimate_errors()
        protein = protein.join(protein_cov)
        return protein

    def __call__(self, timepoints):
        """output: Np x Nt array"""

        with t.no_grad():
            temperature = t.Tensor([self.temperature])
            X = t.Tensor(self.series.cov.X)  # Np x Nr
            k_int = t.Tensor(self.series.cov['k_int'].to_numpy()).unsqueeze(-1)  # Nr x 1
            timepoints = t.Tensor(timepoints).unsqueeze(0)  # 1 x Nt
            inputs = [temperature, X, k_int, timepoints]

            output = self.model(*inputs)
        return output.detach().numpy()
