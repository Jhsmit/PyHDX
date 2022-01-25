from copy import deepcopy

import numpy as np
import pandas as pd
import torch as t
import torch.nn as nn
from scipy import constants, linalg

from pyhdx.fileIO import dataframe_to_file
from pyhdx.models import Protein
from pyhdx.config import cfg


# TORCH_DTYPE = t.double
# TORCH_DEVICE = t.device('cpu')


class DeltaGFit(nn.Module):
    def __init__(self, dG):
        super(DeltaGFit, self).__init__()
        self.register_parameter(name="dG", param=nn.Parameter(dG))

    def forward(self, temperature, X, k_int, timepoints):
        """
        # inputs, list of:
            temperatures: scalar (1,)
            X (N_peptides, N_residues)
            k_int: (N_peptides, 1)

        """

        pfact = t.exp(self.dG / (constants.R * temperature))
        uptake = 1 - t.exp(-t.matmul((k_int / (1 + pfact)), timepoints))
        return t.matmul(X, uptake)


def estimate_errors(hdxm, dG):
    """
    Calculate covariances and uncertainty (perr, experimental)

    Parameters
    ----------
    hdxm : :class:`~pyhdx.models.HDXMeasurement`
    dG : :class:`~numpy.ndarray`
        Array with dG values.

    Returns
    -------

    """
    dtype = t.float64
    joined = pd.concat([dG, hdxm.coverage["exchanges"]], axis=1, keys=["dG", "ex"])
    dG = joined.query("ex==True")["dG"]
    dG_tensor = t.tensor(dG.to_numpy(), dtype=dtype)

    tensors = {
        k: v.cpu() for k, v in hdxm.get_tensors(exchanges=True, dtype=dtype).items()
    }

    def hes_loss(dG_input):
        criterion = t.nn.MSELoss(reduction="sum")
        pfact = t.exp(dG_input.unsqueeze(-1) / (constants.R * tensors["temperature"]))
        uptake = 1 - t.exp(
            -t.matmul((tensors["k_int"] / (1 + pfact)), tensors["timepoints"])
        )
        d_calc = t.matmul(tensors["X"], uptake)

        loss = criterion(d_calc, tensors["d_exp"])
        return loss

    hessian = t.autograd.functional.hessian(hes_loss, dG_tensor)
    hessian_inverse = t.inverse(-hessian)
    covariance = np.sqrt(np.abs(np.diagonal(hessian_inverse)))
    cov_series = pd.Series(covariance, index=dG.index, name="covariance")

    def jac_loss(dG_input):
        pfact = t.exp(dG_input.unsqueeze(-1) / (constants.R * tensors["temperature"]))
        uptake = 1 - t.exp(
            -t.matmul((tensors["k_int"] / (1 + pfact)), tensors["timepoints"])
        )
        d_calc = t.matmul(tensors["X"], uptake)

        residuals = d_calc - tensors["d_exp"]

        return residuals.flatten()

    # https://stackoverflow.com/questions/42388139/how-to-compute-standard-deviation-errors-with-scipy-optimize-least-squares
    jacobian = t.autograd.functional.jacobian(jac_loss, dG_tensor).numpy()

    U, s, Vh = linalg.svd(jacobian, full_matrices=False)
    tol = np.finfo(float).eps * s[0] * max(jacobian.shape)
    w = s > tol
    cov = (Vh[w].T / s[w] ** 2) @ Vh[w]  # robust covariance matrix
    res = jac_loss(dG_tensor).numpy()

    chi2dof = np.sum(res ** 2) / (res.size - dG_tensor.numpy().size)
    cov *= chi2dof
    perr = np.sqrt(np.diag(cov))
    perr_series = pd.Series(perr, index=dG.index, name="perr")

    return cov_series, perr_series


class TorchFitResult(object):
    """
    PyTorch Fit result object.

    Parameters
    ----------

    hdxm_set : :class:`~pyhdx.models.HDXMeasurementSet`
    model
    **metdata

    """

    def __init__(self, hdxm_set, model, losses=None, **metadata):
        self.hdxm_set = hdxm_set
        self.model = model
        self.losses = losses
        self.metadata = metadata
        self.metadata["model_name"] = type(model).__name__
        if losses is not None:
            self.metadata["total_loss"] = self.total_loss
            self.metadata["mse_loss"] = self.mse_loss
            self.metadata["reg_loss"] = self.reg_loss
            self.metadata["regularization_percentage"] = self.regularization_percentage
            self.metadata["epochs_run"] = len(self.losses)

        names = [hdxm.name for hdxm in self.hdxm_set.hdxm_list]

        dfs = [
            self.generate_output(hdxm, self.dG[g_column])
            for hdxm, g_column in zip(self.hdxm_set, self.dG)
        ]
        df = pd.concat(dfs, keys=names, names=["state", "quantity"], axis=1)

        self.output = df

    def get_peptide_mse(self):
        """Get a dataframe with mean squared error per peptide (ie per peptide squared error averaged over time)"""
        squared_errors = self.get_squared_errors()
        dfs = {}
        for mse_sample, hdxm in zip(squared_errors, self.hdxm_set):
            peptide_data = hdxm[0].data
            mse = np.mean(mse_sample, axis=1)
            # Indexing of mse_sum with Np to account for zero-padding
            passthrough_fields = ["start", "end", "sequence"]
            df = peptide_data[passthrough_fields].copy()
            df["peptide_mse"] = mse[: hdxm.Np]
            dfs[hdxm.name] = df

        mse_df = pd.concat(
            dfs.values(), keys=dfs.keys(), names=["state", "quantity"], axis=1
        )

        return mse_df

    def get_residue_mse(self):
        """Get a dataframe with residue mean squared errors

        Errors are from peptide MSE which is subsequently reduced to residue level by weighted averaging

        """

        peptide_mse = self.get_peptide_mse()

        residue_mse_list = []
        for hdxm in self.hdxm_set:
            peptide_mse_values = (
                peptide_mse[hdxm.name, "peptide_mse"].dropna().to_numpy()
            )
            residue_mse_values = hdxm.coverage.Z_norm.T.dot(peptide_mse_values)
            residue_mse = pd.Series(residue_mse_values, index=hdxm.coverage.r_number)
            residue_mse_list.append(residue_mse)

        residue_mse = pd.concat(residue_mse_list, keys=self.hdxm_set.names, axis=1)
        columns = pd.MultiIndex.from_tuples(
            [(name, "residue_mse") for name in self.hdxm_set.names],
            names=["state", "quantity"],
        )
        residue_mse.columns = columns

        return residue_mse

    @property
    def mse_loss(self):
        """obj:`float`: Losses from mean squared error part of Lagrangian"""
        mse_loss = self.losses["mse_loss"].iloc[-1]
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
    def dG(self):
        """output dG as :class:`~pandas.Series` or as :class:`~pandas.DataFrame`

        index is residue numbers
        """

        g_values = self.model.dG.cpu().detach().numpy().squeeze()
        # if g_values.ndim == 1:
        #     dG = pd.Series(g_values, index=self.hdxm_set.coverage.index)
        # else:
        dG = pd.DataFrame(
            g_values.T, index=self.hdxm_set.coverage.index, columns=self.hdxm_set.names
        )

        return dG

    @staticmethod
    def generate_output(hdxm, dG):
        """

        Parameters
        ----------
        hdxm : :class:`~pyhdx.models.HDXMeasurement`
        dG : :class:`~pandas.Series` with r_number as index

        Returns
        -------

        """

        sequence = hdxm.coverage["sequence"].reindex(dG.index)

        out_dict = {"sequence": sequence, "_dG": dG}
        out_dict["dG"] = out_dict["_dG"].copy()
        exchanges = hdxm.coverage["exchanges"].reindex(dG.index, fill_value=False)
        out_dict["dG"][~exchanges] = np.nan
        pfact = np.exp(out_dict["dG"] / (constants.R * hdxm.temperature))
        out_dict["pfact"] = pfact

        k_int = hdxm.coverage["k_int"].reindex(dG.index, fill_value=False)

        k_obs = k_int / (1 + pfact)
        out_dict["k_obs"] = k_obs

        covariance, perr = estimate_errors(hdxm, dG)

        index = pd.Index(hdxm.coverage.r_number, name="r_number")
        df = pd.DataFrame(out_dict, index=dG.index)
        df = df.join(covariance)

        return df

    def to_file(
        self,
        file_path,
        include_version=True,
        include_metadata=True,
        fmt="csv",
        **kwargs,
    ):
        """save only output to file"""
        metadata = self.metadata if include_metadata else include_metadata
        dataframe_to_file(
            file_path,
            self.output,
            include_version=include_version,
            include_metadata=metadata,
            fmt=fmt,
            **kwargs,
        )

    def get_squared_errors(self) -> np.ndarray:
        """np.ndarray: Returns the squared error per peptide per timepoint. Output shape is Ns x Np x Nt"""

        d_calc = self(self.hdxm_set.timepoints)
        errors = (d_calc - self.hdxm_set.d_exp) ** 2

        return errors

    def __call__(self, timepoints):
        """timepoints: shape must be Ns x Nt, or Nt and will be reshaped to Ns x 1 x Nt
        output: Ns x Np x Nt array"""
        # todo fix and tests

        timepoints = np.array(timepoints)
        if timepoints.ndim == 1:
            time_reshaped = np.tile(timepoints, (self.hdxm_set.Ns, 1, 1))
        elif timepoints.ndim == 2:
            Ns, Nt = timepoints.shape
            assert (
                Ns == self.hdxm_set.Ns
            ), "First dimension of 'timepoints' must match the number of samples"
            time_reshaped = timepoints.reshape(Ns, 1, Nt)
        elif timepoints.ndim == 3:
            assert (
                timepoints.shape[0] == self.hdxm_set.Ns
            ), "First dimension of 'timepoints' must match the number of samples"
            time_reshaped = timepoints
        else:
            raise ValueError("Invalid timepoints number of dimensions, must be <=3")

        dtype = t.float64
        with t.no_grad():
            tensors = self.hdxm_set.get_tensors()
            inputs = [tensors[key] for key in ["temperature", "X", "k_int"]]

            time_tensor = t.tensor(time_reshaped, dtype=dtype)
            inputs.append(time_tensor)

            output = self.model(*inputs)

        # todo return as dataframe?
        return output.detach().numpy()

    def __len__(self):
        return self.hdxm_set.Ns


# class TorchSingleFitResult(TorchFitResult):
#     def __init__(self, *args, **kwargs):
#         super(TorchSingleFitResult, self).__init__(*args, **kwargs)
#
#         df = self.generate_output(self.hdxm_set, self.dG)
#         self.output = Protein(df)
#
#     def __call__(self, timepoints):
#         """ timepoints: Nt array (will be unsqueezed to 1 x Nt)
#         output: Np x Nt array"""
#         #todo fix and tests
#         dtype = t.float64
#
#         with t.no_grad():
#             tensors = self.hdxm_set.get_tensors()
#             inputs = [tensors[key] for key in ['temperature', 'X', 'k_int']]
#             inputs.append(t.tensor(timepoints, dtype=dtype).unsqueeze(0))
#
#             output = self.model(*inputs)
#         return output.detach().numpy()
#
#     def __len__(self):
#         return 1


# class TorchBatchFitResult(TorchFitResult):
#     def __init__(self, *args, **kwargs):
#         super(TorchBatchFitResult, self).__init__(*args, **kwargs)


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

    def to_dataframe(self, names=None, field="dG"):
        """convert history of `field` into dataframe.
        names must be given for batch fits with length equal to number of states

        """
        entry = next(iter(self.model_history.values()))
        g = entry[field]
        if g.ndim == 3:
            num_states = entry[field].shape[0]  # G shape is Ns x Nr x 1
            if not len(names) == num_states:
                raise ValueError(
                    f"Number of names provided must be equal to number of states ({num_states})"
                )

            dfs = []
            for i in range(num_states):
                df = pd.DataFrame(
                    {
                        k: v[field].numpy()[i].squeeze()
                        for k, v in self.model_history.items()
                    }
                )
                dfs.append(df)
                full_df = pd.concat(dfs, keys=names, axis=1)
        else:
            full_df = pd.DataFrame(
                {k: v[field].numpy().squeeze() for k, v in self.model_history.items()}
            )

        return full_df
