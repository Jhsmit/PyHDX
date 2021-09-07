from copy import deepcopy

import numpy as np
import pandas as pd
import torch as t
import torch.nn as nn
from scipy import constants, linalg

from pyhdx.fileIO import dataframe_to_file
from pyhdx.models import Protein

TORCH_DTYPE = t.double
TORCH_DEVICE = t.device('cpu')

class DeltaGFit(nn.Module):
    def __init__(self, deltaG):
        super(DeltaGFit, self).__init__()
        #self.deltaG = deltaG
        self.register_parameter(name='deltaG', param=nn.Parameter(deltaG))

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


def estimate_errors(hdxm, deltaG):
    """
    Calculate covariances and uncertainty (perr, experimental)

    Parameters
    ----------
    hdxm : :class:`~pyhdx.models.HDXMeasurement`
    deltaG : :class:`~numpy.ndarray`
        Array with deltaG values.

    Returns
    -------

    """
    joined = pd.concat([deltaG, hdxm.coverage['exchanges']], axis=1, keys=['dG', 'ex'])
    dG = joined.query('ex==True')['dG']
    deltaG = t.tensor(dG.to_numpy(), dtype=TORCH_DTYPE)

    tensors = {k: v.cpu() for k, v in hdxm.get_tensors(exchanges=True).items()}

    def hes_loss(deltaG_input):
        criterion = t.nn.MSELoss(reduction='sum')
        pfact = t.exp(deltaG_input.unsqueeze(-1) / (constants.R * tensors['temperature']))
        uptake = 1 - t.exp(-t.matmul((tensors['k_int'] / (1 + pfact)), tensors['timepoints']))
        output = t.matmul(tensors['X'], uptake)

        loss = criterion(output, tensors['uptake'])
        return loss

    hessian = t.autograd.functional.hessian(hes_loss, deltaG)
    hessian_inverse = t.inverse(-hessian)
    covariance = np.sqrt(np.abs(np.diagonal(hessian_inverse)))
    cov_series = pd.Series(covariance, index=dG.index, name='covariance')

    def jac_loss(deltaG_input):
        pfact = t.exp(deltaG_input.unsqueeze(-1) / (constants.R * tensors['temperature']))
        uptake = 1 - t.exp(-t.matmul((tensors['k_int'] / (1 + pfact)), tensors['timepoints']))
        output = t.matmul(tensors['X'], uptake)

        residuals = (output - tensors['uptake'])

        return residuals.flatten()

    # https://stackoverflow.com/questions/42388139/how-to-compute-standard-deviation-errors-with-scipy-optimize-least-squares
    jacobian = t.autograd.functional.jacobian(jac_loss, deltaG).numpy()

    U, s, Vh = linalg.svd(jacobian, full_matrices=False)
    tol = np.finfo(float).eps * s[0] * max(jacobian.shape)
    w = s > tol
    cov = (Vh[w].T / s[w] ** 2) @ Vh[w]  # robust covariance matrix
    res = jac_loss(deltaG).numpy()

    chi2dof = np.sum(res ** 2) / (res.size - deltaG.numpy().size)
    cov *= chi2dof
    perr = np.sqrt(np.diag(cov))
    perr_series = pd.Series(perr, index=dG.index, name='perr')

    return cov_series, perr_series


class TorchFitResult(object):
    """
    PyTorch Fit result object.

    Parameters
    ----------

    data_obj : :class:`~pyhdx.models.HDXMeasurement` or :class:`~pyhdx.models.HDXMeasurementSet`
    model
    **metdata

    """
    def __init__(self, data_obj, model, losses=None, **metadata):
        self.data_obj = data_obj
        self.model = model
        self.losses = losses
        self.metadata = metadata
        self.metadata['model_name'] = type(model).__name__
        if losses is not None:
            self.metadata['total_loss'] = self.total_loss
            self.metadata['mse_loss'] = self.mse_loss
            self.metadata['reg_loss'] = self.reg_loss
            self.metadata['regularization_percentage'] = self.regularization_percentage
            self.metadata['epochs_run'] = len(self.losses)
        self.output = None  # implemented by subclasses

    @property
    def mse_loss(self):
        """obj:`float`: Losses from mean squared error part of Lagrangian"""
        mse_loss = self.losses['mse_loss'].iloc[-1]
        return float(mse_loss)

    @property
    def total_loss(self):
        """obj:`float`: Total loss value of the Lagrangian"""
        total_loss = self.losses.iloc[-1].sum()
        return float(total_loss)

    @property
    def reg_loss(self):
        """:obj:`float`: Losses from regularization part of Lagrangian"""
        return self.total_loss - self.mse_loss

    @property
    def regularization_percentage(self):
        """:obj:`float`: Percentage part of the total loss that is regularization loss"""
        return (self.reg_loss / self.total_loss) * 100

    @property
    def deltaG(self):
        """output deltaG as :class:`~pandas.Series` or as :class:`~pandas.DataFrame`

        index is residue numbers
        """

        g_values = self.model.deltaG.cpu().detach().numpy().squeeze()
        if g_values.ndim == 1:
            deltaG = pd.Series(g_values, index=self.data_obj.coverage.index)
        else:
            deltaG = pd.DataFrame(g_values.T, index=self.data_obj.coverage.index, columns=self.data_obj.names)

        return deltaG

    @staticmethod
    def generate_output(hdxm, deltaG):
        """

        Parameters
        ----------
        hdxm : :class:`~pyhdx.models.HDXMeasurement`
        deltaG : :class:`~pandas.Series` with r_number as index

        Returns
        -------

        """
        out_dict = {
                    'sequence': hdxm.coverage['sequence'].to_numpy(),
                    '_deltaG': deltaG}
        out_dict['deltaG'] = out_dict['_deltaG'].copy()
        exchanges = hdxm.coverage['exchanges'].reindex(deltaG.index, fill_value=False)
        out_dict['deltaG'][~exchanges] = np.nan
        pfact = np.exp(out_dict['deltaG'] / (constants.R * hdxm.temperature))
        out_dict['pfact'] = pfact

        k_int = hdxm.coverage['k_int'].reindex(deltaG.index, fill_value=False)

        k_obs = k_int / (1 + pfact)
        out_dict['k_obs'] = k_obs

        covariance, perr = estimate_errors(hdxm, deltaG)

        index = pd.Index(hdxm.coverage.r_number, name='r_number')
        df = pd.DataFrame(out_dict, index=index)
        df = df.join(covariance)
        # df = df.join(perr)

        return df

    def to_file(self, file_path, include_version=True, include_metadata=True, fmt='csv', **kwargs):
        metadata = self.metadata if include_metadata else include_metadata
        dataframe_to_file(file_path, self.output.df, include_version=include_version, include_metadata=metadata,
                          fmt=fmt, **kwargs)

    def get_mse(self):
        """np.ndarray: Returns the mean squared error per peptide per timepoint. Output shape is Np x Nt"""

        d_calc = self(self.data_obj.timepoints)
        mse = (d_calc - self.data_obj.d_exp) ** 2

        return mse


class TorchSingleFitResult(TorchFitResult):
    def __init__(self, *args, **kwargs):
        super(TorchSingleFitResult, self).__init__(*args, **kwargs)

        df = self.generate_output(self.data_obj, self.deltaG)
        self.output = Protein(df)

    def __call__(self, timepoints):
        """ timepoints: Nt array (will be unsqueezed to 1 x Nt)
        output: Np x Nt array"""
        #todo fix and tests
        dtype = t.float64

        with t.no_grad():
            tensors = self.data_obj.get_tensors()
            inputs = [tensors[key] for key in ['temperature', 'X', 'k_int']]
            inputs.append(t.tensor(timepoints, dtype=dtype).unsqueeze(0))

            output = self.model(*inputs)
        return output.detach().numpy()


class TorchBatchFitResult(TorchFitResult):
    def __init__(self, *args, **kwargs):
        super(TorchBatchFitResult, self).__init__(*args, **kwargs)
        names = [hdxm.name for hdxm in self.data_obj.hdxm_list]

        dfs = [self.generate_output(hdxm, self.deltaG[g_column]) for hdxm, g_column in zip(self.data_obj, self.deltaG)]
        df = pd.concat(dfs, keys=names, axis=1)

        self.output = Protein(df)

    def __call__(self, timepoints):
        """timepoints: must be Ns x Nt, will be reshaped to Ns x 1 x Nt
        output: Ns x Np x Nt array"""
        #todo fix and tests
        dtype = t.float64
        assert timepoints.shape[0] == self.data_obj.Ns, 'Invalid shape of timepoints'
        with t.no_grad():
            tensors = self.data_obj.get_tensors()
            inputs = [tensors[key] for key in ['temperature', 'X', 'k_int']]

            time_tensor = t.tensor(timepoints.reshape(self.data_obj.Ns, 1, timepoints.shape[1]), dtype=dtype)
            inputs.append(time_tensor)

            output = self.model(*inputs)
        return output.detach().numpy()


class Callback(object):

    def __call__(self, epoch, model, optimizer):
        pass


class CheckPoint(Callback):

    def __init__(self, epoch_step=1000):
        self.epoch_step = epoch_step
        self.model_history = {}

    def __call__(self, epoch, model, optimizer):
        if epoch % self.epoch_step == 0:
            self.model_history[epoch] = deepcopy(model.state_dict())

    def to_dataframe(self, names=None, field='deltaG'):
        """convert history of `field` into dataframe.
         names must be given for batch fits with length equal to number of states

         """
        entry = next(iter(self.model_history.values()))
        g = entry[field]
        if g.ndim == 3:
            num_states = entry[field].shape[0]  # G shape is Ns x Nr x 1
            if not len(names) == num_states:
                raise ValueError(f"Number of names provided must be equal to number of states ({num_states})")

            dfs = []
            for i in range(num_states):
                df = pd.DataFrame({k: v[field].numpy()[i].squeeze() for k, v in self.model_history.items()})
                dfs.append(df)
                full_df = pd.concat(dfs, keys=names, axis=1)
        else:
            full_df = pd.DataFrame({k: v[field].numpy().squeeze() for k, v in self.model_history.items()})

        return full_df