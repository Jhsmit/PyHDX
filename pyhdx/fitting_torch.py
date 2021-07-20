import torch.nn as nn
import torch.nn.functional as F
from torch.optim import SGD
import torch as t
from scipy import constants
import numpy as np
import pandas as pd
from pyhdx.models import Protein
from pyhdx.fileIO import dataframe_to_file


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
    Calculate covariances

    Parameters
    ----------
    hdxm: HDXMeasurement
    deltaG: numpy array
        Array with deltaG values.

    Returns
    -------

    """
    joined = pd.concat([deltaG, hdxm.coverage['exchanges']], axis=1, keys=['dG', 'ex'])
    dG = joined.query('ex==True')['dG']
    deltaG = t.tensor(dG.to_numpy(), dtype=t.float64)

    tensors = hdxm.get_tensors(exchanges=True)

    def calc_loss(deltaG_input):
        criterion = t.nn.MSELoss(reduction='sum')
        pfact = t.exp(deltaG_input.unsqueeze(-1) / (constants.R * tensors['temperature']))
        uptake = 1 - t.exp(-t.matmul((tensors['k_int'] / (1 + pfact)), tensors['timepoints']))
        output = t.matmul(tensors['X'], uptake)

        loss = criterion(output, tensors['uptake'])
        return loss

    hessian = t.autograd.functional.hessian(calc_loss, deltaG)
    hessian_inverse = t.inverse(-hessian)
    covariance = np.sqrt(np.abs(np.diagonal(hessian_inverse)))

    return pd.Series(covariance, index=dG.index, name='covariance')


class TorchFitResult(object):
    """
    PyTorch Fit result object.

    Parameters
    ----------

    data_obj: HDXMeasurement or HDXMeasurementSet
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
        total_loss = self.losses['total_loss'].iloc[-1]
        return float(total_loss)

    @property
    def reg_loss(self):
        """obj:`float`: Losses from regularization part of Lagrangian"""
        return self.total_loss - self.mse_loss

    @property
    def regularization_percentage(self):
        """obj:`float`: Percentage part of the total loss that is regularization loss"""
        return (self.reg_loss / self.total_loss) * 100

    @property
    def deltaG(self):
        """output deltaG as pandas series or as pandas dataframe

        index is residue numbers
        """

        g_values = self.model.deltaG.detach().numpy().squeeze()
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
        hdxm HDXMeasreuement
        deltaG pandas series with r_number as index

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

        covariance = estimate_errors(hdxm, deltaG)

        index = pd.Index(hdxm.coverage.r_number, name='r_number')
        df = pd.DataFrame(out_dict, index=index)
        df = df.join(covariance)

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

        with t.no_grad():
            tensors = self.data_obj.get_tensors()
            inputs = [tensors[key] for key in ['temperature', 'X', 'k_int']]

            time_tensor = t.tensor(timepoints.reshape(self.data_obj.Ns, 1, self.data_obj.Nt), dtype=dtype)
            inputs.append(time_tensor)

            output = self.model(*inputs)
        return output.detach().numpy()



if __name__ == '__main__':
    print('hoi')
    import numpy as np
    import pandas as pd
    data = np.random.rand(100)
    pd_s = pd.Series(data)
    g = t.nn.Parameter(t.tensor(pd_s).unsqueeze(-1))
    print(g.shape)
    module = DeltaGFit(g)
    print(module)
    print(module.__class__)
    print(type(module).__name__)

    g_new = vars()['DeltaGFit'](g)
    print(g_new)