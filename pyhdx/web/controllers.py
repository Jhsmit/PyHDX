import operator
import urllib.request
import zipfile
from collections import namedtuple
from io import StringIO, BytesIO
from pathlib import Path

import colorcet
import dask
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import panel as pn
import param
from numpy.lib.recfunctions import append_fields
from skimage.filters import threshold_multiotsu

from pyhdx import VERSION_STRING
from pyhdx.fileIO import read_dynamx, csv_to_protein, csv_to_dataframe, dataframe_to_stringio
from pyhdx.fitting import fit_rates_weighted_average, fit_rates_half_time_interpolate, get_bounds, fit_gibbs_global, \
    fit_gibbs_global_batch, PATIENCE, STOP_LOSS, EPOCHS, R1, R2, optimizer_defaults
from pyhdx.models import PeptideMasterTable, HDXMeasurement, Protein, array_intersection
from pyhdx.web.base import ControlPanel, DEFAULT_COLORS, DEFAULT_CLASS_COLORS
from pyhdx.web.sources import DataSource, DataFrameSource
from pyhdx.web.transforms import ApplyCmapTransform
from pyhdx.web.widgets import ASyncProgressBar
from pyhdx.support import rgb_to_hex, hex_to_rgba, series_to_pymol

HalfLifeFitResult = namedtuple('HalfLifeFitResult', ['output'])


class MappingFileInputControl(ControlPanel):
    """
    This controller allows users to upload *.txt files where quantities (protection factors, Gibbs free energy, etc) are
    mapped to a linear sequence. The data is then used further downstream to generate binary comparisons between datasets.

    The column should be tab separated with on the last header line (starts with '#') the names of the columns. Columns
    should be tab-delimited.
    """
    header = 'File Input'

    input_file = param.Parameter(default=None, doc='Input file to add to available datasets')
    dataset_name = param.String(doc='Name for the dataset to add. Defaults to filename')
    offset = param.Integer(default=0, doc="Offset to add to the file's r_number column")
    add_dataset = param.Action(lambda self: self._action_add_dataset(),
                               doc='Add the dataset to available datasets')
    datasets_list = param.ListSelector(doc='Current datasets', label='Datasets')
    remove_dataset = param.Action(lambda self: self._action_remove_dataset(),
                                  doc='Remove selected datasets')

    def __init__(self, parent, **params):
        super(MappingFileInputControl, self).__init__(parent, **params)
        self.parent.param.watch(self._datasets_updated, ['datasets'])

    def make_dict(self):
        return self.generate_widgets(input_file=pn.widgets.FileInput)

    @param.depends('input_file', watch=True)
    def _input_file_updated(self):
        self.dataset_name = self.dataset_name or Path(self.widget_dict['input_file'].filename).stem

    @property
    def protein(self):
        """The protein object from the currently selected file in the file widget"""

        try:
            sio = StringIO(self.input_file.decode())
        except UnicodeDecodeError:
            self.parent.logger.info('Invalid file type, supplied file is not a text file')
            return None
        try:
            sio.seek(0)
            protein = txt_to_protein(sio)
        except KeyError:
            sio.seek(0)
            protein = csv_to_protein(sio)
        return protein

    def _add_dataset(self):
        self.parent.datasets[self.dataset_name] = self.protein

    #todo refactor dataset to protein_something
    def _action_add_dataset(self):
        if self.dataset_name in self.parent.datasets.keys():
            self.parent.logger.info(f'Dataset {self.dataset_name} already added')
        elif not self.dataset_name:
            self.parent.logger.info('The added comparison needs to have a name')
        elif not self.input_file:
            self.parent.logger.info('Empty or no file selected')
        elif self.protein is not None:
            self._add_dataset()
            self.parent.param.trigger('datasets')

        self.widget_dict['input_file'].filename = ''
        self.widget_dict['input_file'].value = b''

        self.dataset_name = ''

    def _action_remove_dataset(self):
        if self.datasets_list is not None:
            for dataset_name in self.datasets_list:
                self.parent.datasets.pop(dataset_name)
            self.parent.param.trigger('datasets')

    def _datasets_updated(self, events):
        self.param['datasets_list'].objects = list(self.parent.datasets.keys())


import itertools
cmap_cycle = itertools.cycle(['gray','PiYG', 'jet'])


class CSVFileInputControl(ControlPanel):
    input_file = param.Parameter()
    load_file = param.Action(lambda self: self._action_load())
    temp_new_data = param.Action(lambda self: self._action_new_data())
    temp_new_cmap = param.Action(lambda self: self._action_new_cmap())

    temp_update_filter = param.Action(lambda self: self._action_exposure())
    temp_cmap_rect = param.Action(lambda self: self._action_cmap_rect())

    #cmap_obj = param.ObjectSelector(default='viridis', objects=['viridis', 'plasma', 'magma'])


    def make_dict(self):
        return self.generate_widgets(input_file=pn.widgets.FileInput(accept='.csv,.txt'))

    def _action_load(self):
        sio = StringIO(self.input_file.decode('UTF-8'))
        df = csv_to_dataframe(sio)
        source = DataFrameSource(df=df)

    def _action_new_data(self):

        source = self.parent.sources['torch_fit']
        table = source.get('torch_fit')

        size = len(table)

        new_data = 40e3*np.random.rand(size)

        table['deltaG'] = new_data

        self.parent.update()

    def _action_new_cmap(self):
        cmap_name = np.random.choice(['viridis', 'inferno', 'plasma'])
        cmap = mpl.cm.get_cmap(cmap_name)

        transform = self.parent.transforms['cmap']
        transform.cmap = cmap

        self.parent.update()

    def _action_exposure(self):
        filter = self.parent.filters['exposure']
        filter.widget.value = 0.

        self.parent.update()

    def _action_cmap_rect(self):
        new_cmap = next(cmap_cycle)

        rect_view = self.parent.figure_panels['rect_plot']
        rect_view.opts['cmap'] = new_cmap

        self.parent.update()

        item = self.parent.rows['rect_plot'][0]
        #item.param.trigger('object')


class TestFileInputControl(ControlPanel):
    input_file = param.Parameter()
    load_file = param.Action(lambda self: self._action_load())


    _layout = {
        'self': None,
        'filters.exposure_slider': None
    }

    def __init__(self, parent, **params):
        super().__init__(parent, **params)
        # todo property and list of tuples
        self._layout = {
            'self': None,
            'filters.exposure_slider': None
        }

        self.update_box()

    def make_dict(self):
        return self.generate_widgets(input_file=pn.widgets.FileInput(accept='.csv,.txt'))

    def _action_load(self):
        sio = StringIO(self.input_file.decode('UTF-8'))
        df = csv_to_dataframe(sio)
        source = DataFrameSource(df=df)


class PeptideFileInputControl(ControlPanel):
    """
    This controller allows users to input .csv file (Currently only DynamX format) of 'state' peptide uptake data.
    Users can then choose how to correct for back-exchange and which 'state' and exposure times should be used for
    analysis.

    """
    header = 'Peptide Input'

    input_files = param.List()

    be_mode = param.Selector(doc='Select method of back exchange correction', label='Back exchange correction method', objects=['FD Sample', 'Flat percentage'])
    fd_state = param.Selector(doc='State used to normalize uptake', label='FD State')
    fd_exposure = param.Selector(doc='Exposure used to normalize uptake', label='FD Exposure')
    exp_state = param.Selector(doc='State for selected experiment', label='Experiment State')
    exp_exposures = param.ListSelector(default=[], objects=[''], label='Experiment Exposures'
                                       , doc='Selected exposure time to use')

    be_percent = param.Number(28., bounds=(0, 100), doc='Global percentage of back-exchange',
                              label='Back exchange percentage')

    drop_first = param.Integer(1, bounds=(0, None), doc='Select the number of N-terminal residues to ignore.')
    ignore_prolines = param.Boolean(True, constant=True, doc='Prolines are ignored as they do not exchange D.')
    d_percentage = param.Number(95., bounds=(0, 100), doc='Percentage of deuterium in the labelling buffer',
                                label='Deuterium percentage')
    #fd_percentage = param.Number(95., bounds=(0, 100), doc='Percentage of deuterium in the FD control sample buffer',
    #                             label='FD Deuterium percentage')
    temperature = param.Number(293.15, bounds=(0, 373.15), doc='Temperature of the D-labelling reaction',
                               label='Temperature (K)')
    pH = param.Number(7.5, doc='pH of the D-labelling reaction, as read from pH meter',
                      label='pH read')
    #load_button = param.Action(lambda self: self._action_load(), doc='Load the selected files', label='Load Files')

    n_term = param.Integer(1, doc='Index of the n terminal residue in the protein. Can be set to negative values to '
                                  'accommodate for purification tags. Used in the determination of intrinsic rate of exchange')
    c_term = param.Integer(0, bounds=(0, None),
                           doc='Index of the c terminal residue in the protein. Used for generating pymol export script'
                               'and determination of intrinsic rate of exchange for the C-terminal residue')
    sequence = param.String('', doc='Optional FASTA protein sequence')
    dataset_name = param.String()
    add_dataset_button = param.Action(lambda self: self._action_add_dataset(), label='Add dataset',
                                doc='Parse selected peptides for further analysis and apply back-exchange correction')
    dataset_list = param.ListSelector(label='Datasets', doc='Lists available datasets')

    def __init__(self, parent, **params):
        super(PeptideFileInputControl, self).__init__(parent, **params)
        self.parent.param.watch(self._datasets_updated, ['data_objects'])

        excluded = ['be_percent']
        self.own_widget_names = [name for name in self.widgets.keys() if name not in excluded]
        self.update_box()

        self._array = None  # Numpy array with raw input data

    @property
    def _layout(self):
        return [('self', self.own_widget_names)]

    def make_dict(self):
        text_area = pn.widgets.TextAreaInput(name='Sequence (optional)', placeholder='Enter sequence in FASTA format', max_length=10000,
                                             width=300, height=100, height_policy='fixed', width_policy='fixed')
        return self.generate_widgets(
            input_files=pn.widgets.FileInput(multiple=True, name='Input files'),
            temperature=pn.widgets.FloatInput,
            #be_mode=pn.widgets.RadioButtonGroup,
            be_percent=pn.widgets.FloatInput,
            d_percentage=pn.widgets.FloatInput,
            #fd_percentage=pn.widgets.FloatInput,
            sequence=text_area)

    def make_list(self):
        excluded = ['be_percent']
        widget_list = [widget for name, widget, in self.widget_dict.items() if name not in excluded]

        return widget_list

    @param.depends('be_mode', watch=True)
    def _update_be_mode(self):
        # todo @tejas: Add test
        if self.be_mode == 'FD Sample':
            excluded = ['be_percent']
        elif self.be_mode == 'Flat percentage':
            excluded = ['fd_state', 'fd_exposure']

        self.own_widget_names = [name for name in self.widgets.keys() if name not in excluded]
        #self._layout = {'self': widgets}
        self.update_box()

    @param.depends('input_files', watch=True)
    def _read_files(self):
        """"""
        if self.input_files:
            combined_array = read_dynamx(*[StringIO(byte_content.decode('UTF-8')) for byte_content in self.input_files])
            self._array = combined_array

            self.parent.logger.info(
                f'Loaded {len(self.input_files)} file{"s" if len(self.input_files) > 1 else ""} with a total '
                f'of {len(self._array)} peptides')

        else:
            self._array = None

        self._update_fd_state()
        self._update_fd_exposure()
        self._update_exp_state()
        self._update_exp_exposure()

    def _update_fd_state(self):
        if self._array is not None:
            states = list(np.unique(self._array['state']))
            self.param['fd_state'].objects = states
            self.fd_state = states[0]
        else:
            self.param['fd_state'].objects = []

    @param.depends('fd_state', watch=True)
    def _update_fd_exposure(self):
        if self._array is not None:
            fd_entries = self._array[self._array['state'] == self.fd_state]
            exposures = list(np.unique(fd_entries['exposure']))
        else:
            exposures = []
        self.param['fd_exposure'].objects = exposures
        if exposures:
            self.fd_exposure = exposures[0]

    @param.depends('fd_state', 'fd_exposure', watch=True)
    def _update_exp_state(self):
        if self._array is not None:
            # Booleans of data entries which are in the selected control
            control_bools = np.logical_and(self._array['state'] == self.fd_state, self._array['exposure'] == self.fd_exposure)

            control_data = self._array[control_bools]
            other_data = self._array[~control_bools]

            intersection = array_intersection([control_data, other_data], fields=['start', 'end'])  # sequence?
            states = list(np.unique(intersection[1]['state']))
        else:
            states = []

        self.param['exp_state'].objects = states
        if states:
            self.exp_state = states[0] if not self.exp_state else self.exp_state

    @param.depends('exp_state', watch=True)
    def _update_exp_exposure(self):
        if self._array is not None:
            exp_entries = self._array[self._array['state'] == self.exp_state]
            exposures = list(np.unique(exp_entries['exposure']))
            exposures.sort()
        else:
            exposures = []

        self.param['exp_exposures'].objects = exposures
        self.exp_exposures = exposures

        if not self.dataset_name or self.dataset_name in self.param['exp_state'].objects:
            self.dataset_name = self.exp_state

        if not self.c_term and exposures:
            self.c_term = int(np.max(exp_entries['end']))

    def _datasets_updated(self, events):
        # Update datasets widget as datasets on parents change
        objects = list(self.parent.data_objects.keys())
        self.param['dataset_list'].objects = objects

    def _action_add_dataset(self):
        """Apply controls to :class:`~pyhdx.models.PeptideMasterTable` and set :class:`~pyhdx.models.HDXMeasurement`"""

        if self._array is None:
            self.parent.logger.info("No data loaded")
            return
        elif self.dataset_list and self.dataset_name in self.dataset_list:
            self.parent.logger.info(f"Dataset name {self.dataset_name} already in use")
            return

        peptides = PeptideMasterTable(self._array, d_percentage=self.d_percentage,
                                      drop_first=self.drop_first, ignore_prolines=self.ignore_prolines)
        if self.be_mode == 'FD Sample':
            control_0 = None # = (self.zero_state, self.zero_exposure) if self.zero_state != 'None' else None
            peptides.set_control((self.fd_state, self.fd_exposure), control_0=control_0)
        elif self.be_mode == 'Flat percentage':
            # todo @tejas: Add test
            peptides.set_backexchange(self.be_percent)

        data_states = peptides.data[peptides.data['state'] == self.exp_state]
        data = data_states[np.isin(data_states['exposure'], self.exp_exposures)]

        #todo temperature ph kwarg for series
        hdxm = HDXMeasurement(data, c_term=self.c_term, n_term=self.n_term, sequence=self.sequence,
                              name=self.dataset_name, temperature=self.temperature, pH=self.pH)

        self.parent.data_objects[self.dataset_name] = hdxm
        self.parent.param.trigger('data_objects')  # Trigger update

        df = pd.DataFrame(hdxm.full_data)
        df['start_end'] = [str(s) + '_' + str(e) for s, e in zip(df['start'], df['end'])]
        df['id'] = df.index % hdxm.Np
        target_source = self.parent.sources['dataframe']
        target_source.add_df(df, 'peptides', self.dataset_name)

        index = pd.Index(hdxm.coverage.r_number, name='r_number')
        df = pd.DataFrame(hdxm.rfu_residues, index=index, columns=hdxm.timepoints)
        target_source = self.parent.sources['dataframe']
        target_source.add_df(df, 'rfu', self.dataset_name)

        self.parent.logger.info(f'Loaded dataset {self.dataset_name} with experiment state {self.exp_state} '
                                f'({len(hdxm)} timepoints, {len(hdxm.coverage)} peptides each)')
        self.parent.logger.info(f'Average coverage: {hdxm.coverage.percent_coverage:.3}%, '
                                f'Redundancy: {hdxm.coverage.redundancy:.2}')

    def _action_remove_datasets(self):
        raise NotImplementedError('Removing datasets not implemented')
        for name in self.dataset_list:
            self.parent.datasets.pop(name)

        self.parent.param.trigger('datasets')  # Manual trigger as key assignment does not trigger the param


# todo class DataManagerControl()


class CoverageControl(ControlPanel):
    header = 'Coverage'

    #temp_new_data = param.Action(lambda self: self._action_new_data())

    def __init__(self, parent, **params):
        super().__init__(parent, **params)

        self.update_box()

    @property
    def _layout(self):
        return [
            # ('filters.coverage_state_name', None),
            # ('filters.coverage_exposure', None),
            ('opts.cmap', None),
            #('self', None)
        ]


class InitialGuessControl(ControlPanel):
    """
    This controller allows users to derive initial guesses for D-exchange rate from peptide uptake data.
    """

    #todo remove lambda symbol although its really really funny
    header = 'Initial Guesses'
    fitting_model = param.Selector(default='Half-life (λ)', objects=['Half-life (λ)', 'Association'],
                                   doc='Choose method for determining initial guesses.')
    dataset = param.Selector(default='', doc='Dataset to apply bounds to', label='Dataset (for bounds)')
    global_bounds = param.Boolean(default=False, doc='Set bounds globally across all datasets')
    lower_bound = param.Number(0., doc='Lower bound for association model fitting')
    upper_bound = param.Number(0., doc='Upper bound for association model fitting')
    guess_name = param.String(default='Guess_1', doc='Name for the initial guesses')
    do_fit1 = param.Action(lambda self: self._action_fit(), label='Calculate Guesses', doc='Start initial guess fitting',
                           constant=True)

    bounds = param.Dict({}, doc='Dictionary which stores rate fitting bounds', precedence=-1)

    def __init__(self, parent, **params):
        self.pbar1 = ASyncProgressBar()  #tqdm? https://github.com/holoviz/panel/pull/2079
        self.pbar2 = ASyncProgressBar()
        super(InitialGuessControl, self).__init__(parent, **params)
        self.parent.param.watch(self._parent_datasets_updated, ['data_objects'])  #todo refactor

        excluded = ['lower_bound', 'upper_bound', 'global_bounds', 'dataset']
        self.own_widget_names = [name for name in self.widgets.keys() if name not in excluded]
        self.update_box()

        self._guess_names = {}

    @property
    def _layout(self):
        return [
            ('self', self.own_widget_names),
            # ('filters.select_index_rates_lv1', None),
            # ('filters.select_index_rates_lv2', None),
                        ]

    def make_dict(self):
        widgets = self.generate_widgets(lower_bound=pn.widgets.FloatInput, upper_bound=pn.widgets.FloatInput)
        widgets.update(pbar1=self.pbar1.view, pbar2=self.pbar2.view)

        return widgets

    @param.depends('fitting_model', watch=True)
    def _fitting_model_updated(self):
        if self.fitting_model == 'Half-life (λ)':
            excluded = ['dataset', 'lower_bound', 'upper_bound', 'global_bounds']

        elif self.fitting_model in ['Association', 'Dissociation']:
            excluded = []

        self.own_widget_names = [name for name in self.widgets.keys() if name not in excluded]
        self.update_box()

    @param.depends('global_bounds', watch=True)
    def _global_bounds_updated(self):
        if self.global_bounds:
            self.param['dataset'].constant = True
        else:
            self.param['dataset'].constant = False

    @param.depends('dataset', watch=True)
    def _dataset_updated(self):
        lower, upper = self.bounds[self.dataset]
        self.lower_bound = lower
        self.upper_bound = upper

    @param.depends('lower_bound', 'upper_bound', watch=True)
    def _bounds_updated(self):
        # if self.global_bounds:
        #     for k in self.bounds.keys():
        #         self.bounds[k] = (self.lower_bound, self.upper_bound)
        if not self.global_bounds:
            self.bounds[self.dataset] = (self.lower_bound, self.upper_bound)

    def _parent_datasets_updated(self, events):
        if len(self.parent.data_objects) > 0:
            self.param['do_fit1'].constant = False

        # keys to remove:
        for k in self.bounds.keys() - self.parent.data_objects.keys():
            self.bounds.pop(k)
        # keys to add:
        for k in self.parent.data_objects.keys() - self.bounds.keys():
            self.bounds[k] = get_bounds(self.parent.data_objects[k].timepoints)

        options = list(self.parent.data_objects.keys())
        self.param['dataset'].objects = options
        if not self.dataset:
            self.dataset = options[0]

    def add_fit_result(self, future):
        name = self._guess_names.pop(future.key)

        results = future.result()
        dfs = [result.output.df for result in results]
        combined_results = pd.concat(dfs, axis=1,
                                     keys=list(self.parent.data_objects.keys()),
                                     names=['state_name', 'quantity'])

        self.sources['dataframe'].add_df(combined_results, 'rates', name)
        self.parent.fit_results[name] = {k: v for k, v in zip(self.parent.data_objects.keys(), results)}
        self.parent.param.trigger('data_objects')  # Informs other fittings that initial guesses are now available
        self.param['do_fit1'].constant = False

    def _action_fit(self):
        if len(self.parent.data_objects) == 0:
            self.parent.logger.info('No datasets loaded')
            return

        if self.guess_name in itertools.chain(self.parent.fit_results.keys(), self._guess_names.values()):
            self.parent.logger.info(f"Guess with name {self.guess_name} already in use")
            return

        self.parent.logger.debug('Start initial guess fit')
        self.param['do_fit1'].constant = True

        num_samples = len(self.parent.data_objects)
        if self.fitting_model.lower() in ['association', 'dissociation']:
            if self.global_bounds:
                bounds = [(self.lower_bound, self.upper_bound)]*num_samples
            else:
                bounds = self.bounds.values()

            futures = self.parent.client.map(fit_rates_weighted_average,
                                             self.parent.data_objects.values(), bounds, client='worker_client')
        elif self.fitting_model == 'Half-life (λ)':   # this is practically instantaneous and does not require dask
            futures = self.parent.client.map(fit_rates_half_time_interpolate, self.parent.data_objects.values())

        dask_future = self.parent.client.submit(lambda args: args, futures)
        self._guess_names[dask_future.key] = self.guess_name

        self.parent.future_queue.append((dask_future, self.add_fit_result))


class FitControl(ControlPanel):
    """
    This controller allows users to execute PyTorch fitting of the global data set.

    Currently, repeated fitting overrides the old result.
    """

    header = 'Fitting'

    initial_guess = param.Selector(doc='Name of dataset to use for initial guesses.')

    fit_mode = param.Selector(default='Batch', objects=['Batch', 'Single'])

    stop_loss = param.Number(STOP_LOSS, bounds=(0, None),
                             doc='Threshold loss difference below which to stop fitting.')
    stop_patience = param.Integer(PATIENCE, bounds=(1, None),
                                  doc='Number of epochs where stop loss should be satisfied before stopping.')
    learning_rate = param.Number(optimizer_defaults['SGD']['lr'], bounds=(0, None),
                                 doc='Learning rate parameter for optimization.')
    momentum = param.Number(optimizer_defaults['SGD']['momentum'], bounds=(0, None),
                            doc='Stochastic Gradient Descent momentum')
    nesterov = param.Boolean(optimizer_defaults['SGD']['nesterov'],
                             doc='Use Nesterov type of momentum for SGD')
    epochs = param.Integer(EPOCHS, bounds=(1, None),
                           doc='Maximum number of epochs (iterations.')
    r1 = param.Number(R1, bounds=(0, None), label='Regularizer 1 (peptide axis)',
                      doc='Value of the regularizer along residue axis.')

    r2 = param.Number(R2, bounds=(0, None), label='Regularizer 2 (sample axis)',
                      doc='Value of the regularizer along sample axis.', constant=True)

    fit_name = param.String("Gibbs_fit_1", doc="Name for for the fit result")

    do_fit = param.Action(lambda self: self._action_fit(), constant=True, label='Do Fitting',
                          doc='Start global fitting')

    def __init__(self, parent, **params):
        self.pbar1 = ASyncProgressBar() #tqdm?
        super(FitControl, self).__init__(parent, **params)

        source = self.parent.sources['dataframe']
        source.param.watch(self._source_updated, ['updated'])

        self._current_jobs = 0
        self._max_jobs = 2  #todo config
        self._fit_names = {}

    def _source_updated(self, *events):
        table = self.parent.sources['dataframe'].get('rates')

        objects = list(table.columns.levels[0])
        if objects:
            self.param['do_fit'].constant = False

        self._fit_mode_updated()

        self.param['initial_guess'].objects = objects
        if not self.initial_guess and objects:
            self.initial_guess = objects[0]

    @param.depends('fit_mode', watch=True)
    def _fit_mode_updated(self):
        if self.fit_mode == 'Batch' and len(self.parent.data_objects) > 1:
            self.param['r2'].constant = False
        else:
            self.param['r2'].constant = True

    def add_fit_result(self, future):
        #todo perhaps all these dfs should be in the future?
        name = self._fit_names.pop(future.key)
        result = future.result()
        self._current_jobs -= 1

        self.parent.logger.info(f'Finished PyTorch fit: {name}')

        # List of single fit results
        if isinstance(result, list):
            self.parent.fit_results[name] = list(result)
            output_dfs = {fit_result.data_obj.name: fit_result.output.df for fit_result in result}
            df = pd.concat(output_dfs.values(), keys=output_dfs.keys(), axis=1)

            # create mse losses dataframe
            dfs = {}
            for single_result in result:
            # Determine mean squared errors per peptide, summed over timepoints
                mse = single_result.get_mse()
                mse_sum = np.sum(mse, axis=1)
                peptide_data = single_result.data_obj[0].data
                data_dict = {'start': peptide_data['start'], 'end': peptide_data['end'], 'total_mse': mse_sum}
                dfs[single_result.data_obj.name] = pd.DataFrame(data_dict)
            mse_df = pd.concat(dfs.values(), keys=dfs.keys(), axis=1)

        #todo d calc for single fits
        #todo losses for single fits

            # Create d_calc dataframe
            # -----------------------
            # todo needs cleaning up
            state_dfs = {}
            for single_result in result:
                tp_flat = single_result.data_obj.timepoints
                elem = tp_flat[np.nonzero(tp_flat)]

                time_vec = np.logspace(np.log10(elem.min()) - 1, np.log10(elem.max()), num=100, endpoint=True)
                d_calc_state = single_result(time_vec)  #shape Np x Nt
                hdxm = single_result.data_obj

                peptide_dfs = []
                pm_data = hdxm[0].data
                for d_peptide, pm_row in zip(d_calc_state, pm_data):
                    peptide_id = f"{pm_row['start']}_{pm_row['end']}"
                    data_dict = {'timepoints': time_vec, 'd_calc': d_peptide, 'start_end': [peptide_id] * len(time_vec)}
                    peptide_dfs.append(pd.DataFrame(data_dict))
                state_dfs[hdxm.name] = pd.concat(peptide_dfs, axis=0, ignore_index=True)

            d_calc_df = pd.concat(state_dfs.values(), keys=state_dfs.keys(), axis=1)


            # Create losses/epoch dataframe
            # -----------------------------
            losses_dfs = {fit_result.data_obj.name: fit_result.losses for fit_result in result}
            losses_df = pd.concat(losses_dfs.values(), keys=losses_dfs.keys(), axis=1)


        else:  # one batchfit result
            self.parent.fit_results[name] = result  # todo this name can be changed by the time this is executed
            df = result.output.df
            # df.index.name = 'peptide index'

            # Create MSE losses df (per peptide, summed over timepoints)
            # -----------------------
            mse = result.get_mse()
            dfs = {}
            for mse_sample, hdxm in zip(mse, result.data_obj):
                peptide_data = hdxm[0].data
                mse_sum = np.sum(mse_sample, axis=1)
                # Indexing of mse_sum with Np to account for zero-padding
                data_dict = {'start': peptide_data['start'], 'end': peptide_data['end'], 'total_mse': mse_sum[:hdxm.Np]}
                dfs[hdxm.name] = pd.DataFrame(data_dict)

            mse_df = pd.concat(dfs.values(), keys=dfs.keys(), axis=1)

            self.parent.logger.info('Finished PyTorch fit')

            # Create d_calc dataframe
            # -----------------------
            tp_flat = result.data_obj.timepoints.flatten()
            elem = tp_flat[np.nonzero(tp_flat)]

            time_vec = np.logspace(np.log10(elem.min()) - 1, np.log10(elem.max()), num=100, endpoint=True)
            stacked = np.stack([time_vec for i in range(result.data_obj.Ns)])
            d_calc = result(stacked)

            state_dfs = {}
            for hdxm, d_calc_state in zip(result.data_obj, d_calc):
                peptide_dfs = []
                pm_data = hdxm[0].data
                for d_peptide, pm_row in zip(d_calc_state, pm_data):
                    peptide_id = f"{pm_row['start']}_{pm_row['end']}"
                    data_dict = {'timepoints': time_vec, 'd_calc': d_peptide, 'start_end': [peptide_id] * len(time_vec)}
                    peptide_dfs.append(pd.DataFrame(data_dict))
                state_dfs[hdxm.name] = pd.concat(peptide_dfs, axis=0, ignore_index=True)
            d_calc_df = pd.concat(state_dfs.values(), keys=state_dfs.keys(), axis=1)

            # Create losses/epoch dataframe
            # -----------------------------
            losses_df = result.losses.copy()
            losses_df.columns = pd.MultiIndex.from_product(
                [['All states'], losses_df.columns],
                names=['state_name', 'quantity']
            )

            self.parent.logger.info(
                f"Finished fitting in {len(result.losses)} epochs, final mean squared residuals is {result.mse_loss:.2f}")
            self.parent.logger.info(f"Total loss: {result.total_loss:.2f}, regularization loss: {result.reg_loss:.2f} "
                                    f"({result.regularization_percentage:.1f}%)")

        self.parent.sources['dataframe'].add_df(df, 'global_fit', names=[name])
        self.parent.sources['dataframe'].add_df(mse_df, 'peptides_mse', names=[name])
        self.parent.sources['dataframe'].add_df(d_calc_df, 'd_calc', names=[name])
        self.parent.sources['dataframe'].add_df(losses_df, 'losses', names=[name])



        self.parent.param.trigger('fit_results')

    def _action_fit(self):
        if self.fit_name in itertools.chain(self.parent.fit_results.keys(), self._fit_names.values()):
            self.parent.logger.info(f"Fit result with name {self.fit_name} already in use")
            return

        self.parent.logger.info('Started PyTorch fit')

        self._current_jobs += 1
        if self._current_jobs >= self._max_jobs:
            self.widgets['do_fit'].constant = True

        self.parent.logger.info(f'Current number of active jobs: {self._current_jobs}')
        if self.fit_mode == 'Batch':
            hdx_set = self.parent.hdx_set
            rates_df = self.sources['dataframe'].get('rates', fit_ID=self.initial_guess)

            rates_guess = [rates_df[state]['rate'] for state in hdx_set.names]
            gibbs_guess = hdx_set.guess_deltaG(rates_guess)

            dask_future = self.parent.client.submit(fit_gibbs_global_batch, hdx_set, gibbs_guess, **self.fit_kwargs)
        else:
            data_objs = self.parent.data_objects.values()
            rates_df = self.sources['dataframe'].get('rates', fit_ID=self.initial_guess)
            gibbs_guesses = [data_obj.guess_deltaG(rates_df[data_obj.name]['rate']) for data_obj in data_objs]
            futures = self.parent.client.map(fit_gibbs_global, data_objs, gibbs_guesses, **self.fit_kwargs)

            # Combine list of futures into one future object
            # See https://github.com/dask/distributed/pull/560
            dask_future = self.parent.client.submit(lambda args: args, futures)

        self._fit_names[dask_future.key] = self.fit_name
        self.parent.future_queue.append((dask_future, self.add_fit_result))

    @property
    def fit_kwargs(self):
        fit_kwargs = dict(r1=self.r1, lr=self.learning_rate, momentum=self.momentum, nesterov=self.nesterov,
                          epochs=self.epochs, patience=self.stop_patience, stop_loss=self.stop_loss)
        if self.fit_mode == 'Batch':
            fit_kwargs['r2'] = self.r2

        return fit_kwargs


class ClassificationControl(ControlPanel):
    """
    This controller allows users classify 'mapping' datasets and assign them colors.

    Coloring can be either in discrete categories or as a continuous custom color map.
    """

    header = 'Classification'
    # format ['tag1', ('tag2a', 'tag2b') ] = tag1 OR (tag2a AND tag2b)

    # todo unify name for target field (target_data set)
    # When coupling param with the same name together there should be an option to exclude this behaviour
    table = param.Selector(label='Target table')
    # fit_ID = param.Selector()  # generalize selecting widgets based on selected table
    # quantity = param.Selector(label='Quantity')  # this is the lowest-level quantity of the multiindex df (filter??)

    mode = param.Selector(default='Discrete', objects=['Discrete', 'Continuous', 'Color map'],
                          doc='Choose color mode (interpolation between selected colors).')#, 'ColorMap'])
    num_colors = param.Integer(3, bounds=(1, 10), label='Number of colours',
                              doc='Number of classification colors.')
    library = param.Selector(default='matplotlib', objects=['matplotlib', 'colorcet'])
    color_map = param.Selector()
    otsu_thd = param.Action(lambda self: self._action_otsu(), label='Otsu',
                            doc="Automatically perform thresholding based on Otsu's method.")
    linear_thd = param.Action(lambda self: self._action_linear(), label='Linear',
                              doc='Automatically perform thresholding by creating equally spaced sections.')
    log_space = param.Boolean(False,
                              doc='Boolean to set whether to apply colors in log space or not.')
    #apply = param.Action(lambda self: self._action_apply())
    no_coverage = param.Color(default='#8c8c8c', doc='Color to use for regions of no coverage')

    color_set_name = param.String('', doc='Name for the color dataset to add')
    add_colorset = param.Action(lambda self: self._action_add_colorset())

    #show_thds = param.Boolean(True, label='Show Thresholds', doc='Toggle to show/hide threshold lines.')
    values = param.List(default=[], precedence=-1)
    colors = param.List(default=[], precedence=-1)

    def __init__(self, parent, **param):
        super(ClassificationControl, self).__init__(parent, **param)

        # https://discourse.holoviz.org/t/based-on-a-select-widget-update-a-second-select-widget-then-how-to-link-the-latter-to-a-reactive-plot/917/8
        cc_cmaps = sorted(colorcet.cm.keys())
        mpl_cmaps = sorted(set(plt.colormaps()) - set('cet_' + cmap for cmap in cc_cmaps))
        self.cmaps = {'matplotlib': mpl_cmaps, 'colorcet': cc_cmaps}
        self.param['color_map'].objects = mpl_cmaps

        self._update_num_colors()
        self._update_num_values()
        self.excluded = ['library', 'color_map'] # excluded widgets based on choice of `mode`

        views = [view for view in self.views.values() if any(isinstance(trs, ApplyCmapTransform) for trs in view.transforms)]
        options = [view.table for view in views]

        for view in views:
            view.source.param.watch(self._sources_updated, 'updated')

        self.param['table'].objects = options
        if not self.table and options:
            self.table = options[0]

        self._table_updated()  # also updates box
        #self.update_box()

    @property
    def own_widget_names(self):
        """returns a list of names of widgets in self.widgets to be laid out in controller card"""

        # initial_widgets = [name for name in self.widgets.keys() if name not in self.excluded]
        initial_widgets = []
        for name in self.param:
            precedence = self.param[name].precedence
            if (precedence is None or precedence > 0) and name not in self.excluded + ['name']:
                initial_widgets.append(name)
        #l1[1:1] = l2
        select_widgets = [name for name in self.widgets.keys() if name.startswith('select')]
        initial_widgets[1:1] = select_widgets

        #value_widget_names = [f'value_{i}' for i in range(len(self.values))]
        #color_widget_names = [f'color_{i}' for i in range(len(self.colors))]
        widget_names = initial_widgets + [f'value_{i}' for i in range(len(self.values))]
        if self.mode != 'Color map':
            widget_names += [f'color_{i}' for i in range(len(self.colors))]
        return widget_names

    #        return initial_widgets + #list(self.values_widgets.keys()) + list(self.colors_widgets.keys())

    def make_dict(self):
        return self.generate_widgets(num_colors=pn.widgets.IntInput)

    @property
    def _layout(self):
        return [
            ('self', self.own_widget_names),
                ]

    def _sources_updated(self, *events):
        self._table_updated()

    @param.depends('table', watch=True)
    def _table_updated(self):
        df = self.get_data()

        #todo also get schema and check if this table is compatible (ie has numbers, not colors only)
        if df.empty:
            return
        names = df.columns.names

        # Remove old widgets (list comprehension)
        old_widget_names = [key for key in self.widgets.keys() if key.startswith('select')]
        [self.widgets.pop(key) for key in old_widget_names]

        widget_dict = {}
        for i, (name, options) in enumerate(zip(names, df.columns.levels)):
            _opts = ['*'] + list(options) if i != len(names) - 1 else list(options)
            #todo make function to determine defaults
            if i == 0:
                default = _opts[-1]
            else:
                default = 'deltaG' if 'deltaG' in _opts else _opts[0]
            widget = pn.widgets.Select(name=name, options=_opts, value=default)
            widget_dict[f'select_{i}'] = widget

        self.widgets.update(widget_dict)
        self.update_box()

    def get_data(self):
        """object pandas dataframe: returns current multindex dataframe"""
        source = self.sources['dataframe']
        df = source.get(self.table)

        return df

    def get_selected_data(self):
        #todo move method to data source?
        df = self.get_data()
        selected_fields = [widget.value for name, widget in self.widgets.items() if name.startswith('select')]
        bools_list = [df.columns.get_level_values(i) == value for i, value in enumerate(selected_fields) if
                      value != '*']

        if len(bools_list) == 0:
            bools = np.ones(len(df.columns)).astype(bool)
        elif len(bools_list) == 1:
            bools = np.array(bools_list).flatten()
        else:
            bools_array = np.array(bools_list)
            bools = np.product(bools_array, axis=0).astype(bool)

        selected_df = df.iloc[:, bools]

        return selected_df

    def get_values(self):
        """return numpy array with only the values from selected dataframe, nan omitted"""

        array = self.get_selected_data().to_numpy().flatten()
        values = array[~np.isnan(array)]

        return values

    def _action_otsu(self):
        if self.num_colors <= 1:
            return
        values = self.get_values() # todo check for no values
        if not values.size:
            return

        func = np.log if self.log_space else lambda x: x  # this can have NaN when in log space
        thds = threshold_multiotsu(func(values), classes=self.num_colors)
        widgets = [widget for name, widget in self.widgets.items() if name.startswith('value')]
        for thd, widget in zip(thds[::-1], widgets):  # Values from high to low
            widget.start = None
            widget.end = None
            widget.value = np.exp(thd) if self.log_space else thd
        self._update_bounds()

        #self._get_colors()

    def _action_linear(self):
        i = 1 if self.mode == 'Discrete' else 0
        values = self.get_values()
        if not values.size:
            return

        if self.log_space:
            thds = np.logspace(np.log(np.min(values)), np.log(np.max(values)),
                               num=self.num_colors + i, endpoint=True, base=np.e)
        else:
            thds = np.linspace(np.min(values), np.max(values), num=self.num_colors + i, endpoint=True)

        widgets = [widget for name, widget in self.widgets.items() if name.startswith('value')]
        for thd, widget in zip(thds[i:self.num_colors][::-1], widgets):
            # Remove bounds, set values, update bounds
            widget.start = None
            widget.end = None
            widget.value = thd
        self._update_bounds()

    def _action_add_colorset(self):
        if not self.color_set_name:
            self.parent.logger.info('No name given tot the colorset')
            return

        source = self.sources['dataframe']
        if self.color_set_name in source.tables.keys():  #todo update
            self.parent.logger.info(f'Colorset with name {self.color_set_name} already present')
            return

        selected_df = self.get_selected_data()
        cmap, norm = self.get_cmap_and_norm()

        array = cmap(norm(selected_df), bytes=True)
        colors_hex = rgb_to_hex(array.reshape(-1, 4))
        output = colors_hex.reshape(array.shape[:-1])

        output_df = pd.DataFrame(output, index=selected_df.index, columns=selected_df.columns)
        if output_df.index.name == 'r_number':  # The selected dataset is a protein mappable
            c_term = max([data_obj.coverage.protein.c_term for data_obj in self.parent.data_objects.values()])
            n_term = min([data_obj.coverage.protein.n_term for data_obj in self.parent.data_objects.values()])

            new_index = pd.RangeIndex(start=n_term, stop=c_term, name='r_number')
            output_df = output_df.reindex(index=new_index, fill_value=self.no_coverage.upper())
            output_df.rename_axis(columns={'fit_ID': 'color_ID'}, inplace=True)
            output_df.columns = output_df.columns.set_levels([self.color_set_name], level=0)

        source.add_df(output_df, 'colors')

    @param.depends('color_map', 'values', 'colors', watch=True)
    def _action_apply(self):
        cmap, norm = self.get_cmap_and_norm()

        if cmap and norm:
            #this needs to be updated to more generalized
            transform = self.transforms['cmap_transform']
            transform.cmap = cmap
            transform.norm = norm

    def get_cmap_and_norm(self):
        norm_klass = mpl.colors.Normalize if not self.log_space else mpl.colors.LogNorm
        if len(self.values) < 2:
            return None, None

        if self.mode == 'Discrete':
            if len(self.values) != len(self.colors) - 1:
                return None, None
            cmap = mpl.colors.ListedColormap(self.colors)
            norm = mpl.colors.BoundaryNorm(self.values[::-1], self.num_colors, extend='both') #todo refactor values to thd_values
        elif self.mode == 'Continuous':
            norm = norm_klass(vmin=np.min(self.values), vmax=np.max(self.values), clip=True)
            positions = norm(self.values[::-1])
            cmap = mpl.colors.LinearSegmentedColormap.from_list('custom_cmap', list(zip(positions, self.colors)))
        elif self.mode == 'Color map':
            norm = norm_klass(vmin=np.min(self.values), vmax=np.max(self.values), clip=True)
            if self.library == 'matplotlib':
                cmap = mpl.cm.get_cmap(self.color_map)
            elif self.library == 'colorcet':
                cmap = getattr(colorcet, 'm_' + self.color_map)

        cmap.set_bad(self.no_coverage)
        return cmap, norm

    @param.depends('library', watch=True)
    def _update_library(self):
        options = self.cmaps[self.library]
        self.param['color_map'].objects = options

    @param.depends('mode', watch=True)
    def _mode_updated(self):
        if self.mode == 'Discrete':
            self.excluded = ['library', 'color_map']
    #        self.num_colors = max(3, self.num_colors)
    #        self.param['num_colors'].bounds = (3, None)
        elif self.mode == 'Continuous':
            self.excluded = ['library', 'color_map', 'otsu_thd']
      #      self.param['num_colors'].bounds = (2, None)
        elif self.mode == 'Color map':
            self.excluded = ['otsu_thd', 'num_colors']
            self.num_colors = 2

        #todo adjust add/ remove color widgets methods
        self.param.trigger('num_colors')
        self.update_box()

    @param.depends('num_colors', watch=True)
    def _update_num_colors(self):
        while len(self.colors) != self.num_colors:
            if len(self.colors) > self.num_colors:
                self._remove_color()
            elif len(self.colors) < self.num_colors:
                self._add_color()
        self.param.trigger('colors')

    @param.depends('num_colors', watch=True)
    def _update_num_values(self):
        diff = 1 if self.mode == 'Discrete' else 0
        while len(self.values) != self.num_colors - diff:
            if len(self.values) > self.num_colors - diff:
                self._remove_value()
            elif len(self.values) < self.num_colors - diff:
                self._add_value()

        self._update_bounds()
        self.param.trigger('values')
        self.update_box()

    def _add_value(self):
        # value widgets are ordered in decreasing order, ergo next value widget
        # starts with default value of previous value -1
        try:
            first_value = self.values[-1]
        except IndexError:
            first_value = 0

        default = float(first_value - 1)
        self.values.append(default)

        name = f'Threshold {len(self.values)}'
        key = f'value_{len(self.values) - 1}'   # values already populated, first name starts at 1
        widget = pn.widgets.FloatInput(name=name, value=default)
        self.widgets[key] = widget
        widget.param.watch(self._value_event, ['value'])

    def _remove_value(self):
        key = f'value_{len(self.values) - 1}'
        widget = self.widgets.pop(key)
        self.values.pop()

        [widget.param.unwatch(watcher) for watcher in widget.param._watchers]
        del widget

    def _add_color(self):
        try:
            default = DEFAULT_CLASS_COLORS[len(self.colors)]
        except IndexError:
            default = "#"+''.join(np.random.choice(list('0123456789abcdef'), 6))

        self.colors.append(default)

        key = f'color_{len(self.colors) - 1}'
        widget = pn.widgets.ColorPicker(value=default)

        self.widgets[key] = widget

        widget.param.watch(self._color_event, ['value'])

    def _remove_color(self):
        key = f'color_{len(self.colors) - 1}'
        widget = self.widgets.pop(key)
        self.colors.pop()
        [widget.param.unwatch(watcher) for watcher in widget.param._watchers]
        del widget

    def _color_event(self, *events):
        for event in events:
            idx = list(self.widgets.values()).index(event.obj)
            key = list(self.widgets.keys())[idx]
            widget_index = int(key.split('_')[1])
            # idx = list(self.colors_widgets).index(event.obj)
            self.colors[widget_index] = event.new

        self.param.trigger('colors')

        #todo param trigger colors????

    def _value_event(self, *events):
        """triggers when a single value gets changed"""
        for event in events:
            idx = list(self.widgets.values()).index(event.obj)
            key = list(self.widgets.keys())[idx]
            widget_index = int(key.split('_')[1])
            self.values[widget_index] = event.new

        self._update_bounds()
        self.param.trigger('values')

    def _update_bounds(self):
        #for i, widget in enumerate(self.values_widgets.values()):
        for i in range(len(self.values)):
            widget = self.widgets[f'value_{i}']
            if i > 0:
                key = f'value_{i-1}'
                prev_value = float(self.widgets[key].value)
                widget.end = np.nextafter(prev_value, prev_value - 1)
            else:
                widget.end = None

            if i < len(self.values) - 1:
                key = f'value_{i+1}'
                next_value = float(self.widgets[key].value)
                widget.start = np.nextafter(next_value, next_value + 1)
            else:
                widget.start = None


class ProteinControl(ControlPanel):
    header = 'Protein Control'

    input_mode = param.Selector(doc='Method of protein structure input', objects=['PDB File', 'RCSB Download'])
    file_binary = param.Parameter()
    rcsb_id = param.String(doc='RCSB ID of protein to download')
    load_structure = param.Action(lambda self: self._action_load_structure())

    def __init__(self, parent, **params):
        super(ProteinControl, self).__init__(parent, **params)

        excluded = ['rcsb_id']
        self.own_widget_names = [name for name in self.widgets.keys() if name not in excluded]
        self.update_box()

    @property
    def _layout(self):
        return [('self', self.own_widget_names),
                ('filters.ngl_color_id', None),
                ('filters.ngl_state_name', None),
                ]

    def make_dict(self):
        return self.generate_widgets(file_binary=pn.widgets.FileInput(multiple=False, accept='.pdb'))

    @param.depends('input_mode', watch=True)
    def _update_input_mode(self):
        if self.input_mode == 'PDB File':
            excluded = ['rcsb_id']
        elif self.input_mode == 'RCSB Download':
            excluded = ['file_binary']

        self.own_widget_names = [name for name in self.widgets.keys() if name not in excluded]
        self.update_box()

    def _action_load_structure(self):
        view = self.views['protein']
        if self.input_mode == 'PDB File':
            pdb_string = self.file_binary.decode()
            view.ngl_view.pdb_string = pdb_string
        elif self.input_mode == 'RCSB Download':
            if len(self.rcsb_id) != 4:
                self.parent.logger.info(f"Invalid RCSB pdb id: {self.rcsb_id}")
                return

            url = f'http://files.rcsb.org/download/{self.rcsb_id}.pdb'
            with urllib.request.urlopen(url) as response:
                pdb_string = response.read().decode()
                view.ngl_view.pdb_string = pdb_string


class GraphControl(ControlPanel):
    header = 'Graph Control'

    spin = param.Boolean(default=False, doc='Spin the protein object')

    state_name = param.Selector(doc="Name of the currently selected state")
    fit_id = param.Selector(doc="Name of the currently selected fit ID")
    peptide_index = param.Selector(doc="Index of the currently selected peptide")

    def __init__(self, parent, **params):
        super(GraphControl, self).__init__(parent, **params)
        source = self.sources['dataframe']
        source.param.watch(self._source_updated, 'updated')

    def make_dict(self):
        widgets = {
            'general': pn.pane.Markdown('### General'),
            'coverage': pn.pane.Markdown('### Coverage'),
            'peptide': pn.pane.Markdown('### Peptide'),
            'losses': pn.pane.Markdown('### Losses'),
            'debugging': pn.pane.Markdown('### Debugging'),

        }

        return {**widgets, **self.generate_widgets()}

    def _source_updated(self, *events):
        source = self.sources['dataframe']
        table = source.get('global_fit')
        fit_id_options = list(table.columns.get_level_values(0).unique())
        self.param['fit_id'].objects = fit_id_options
        if not self.fit_id and fit_id_options:
            self.fit_id = fit_id_options[0]

        table = source.get('peptides')
        state_name_options = list(table.columns.get_level_values(0).unique())

        self.param['state_name'].objects = state_name_options
        if not self.state_name and state_name_options:
            self.state_name = state_name_options[0]

    @param.depends('state_name', watch=True)
    def _update_state_name(self):
        #https://param.holoviz.org/reference.html#param.parameterized.batch_watch

        dwarfs = ['coverage_state_name', 'coverage_mse_state_name', 'peptide_d_exp_state_name', 'peptide_d_calc_state_name',
                  'deltaG_state_name', 'rates_state_name', 'ngl_state_name']  # there really are 7

        # one filter to rule them all, one filter to find them,
        # one filter to bring them all, and in the darkness bind them;
        # in the Land of Mordor where the shadows lie.
        for dwarf in dwarfs:
            filt = self.filters[dwarf]
            filt.value = self.state_name

        # If current fit result was done as single, also update the state for the losses graph
        losses_filt = self.filters['losses_state_name']
        if self.state_name in losses_filt.param['value'].objects:
            losses_filt.value = self.state_name


        # Update possible choices for peptide selection depending on selected state
        source = self.sources['dataframe']
        table = source.get('peptides')
        unique_vals = table[self.state_name]['start_end'].unique()
        peptide_options = list(range(len(unique_vals)))
        self.param['peptide_index'].objects = peptide_options
        if self.peptide_index is not None and peptide_options:
            self.peptide_index = peptide_options[0]

    @param.depends('fit_id', watch=True)
    def _update_fit_id(self):
        elves = ['coverage_mse_fit_id', 'peptide_d_calc_fit_id', 'deltaG_fit_id', 'losses_fit_id']
        for elf in elves:
            filt = self.filters[elf]
            filt.value = self.fit_id

        # perhaps this is faster?
        # widget = self.widget.clone()
        # self.widget.link(widget, value='value', bidirectional=True)

    @param.depends('peptide_index', watch=True)
    def _update_peptide_index(self):
        hobbits = ['peptide_d_exp_select', 'peptide_d_calc_select']
        for hobbit in hobbits:
            filt = self.filters[hobbit]
            filt.value = self.peptide_index

    @property
    def _layout(self):
        return [
            # ('self', ['coverage']),
            # ('filters.select_index', None),
            # ('filters.exposure_slider', None),
            # ('opts.cmap', None),
            ('self', ['general']),
            ('self', ['fit_id', 'state_name']),
            ('self', ['coverage']),
            ('filters.coverage_exposure', None),
            ('self', ['peptide', 'peptide_index']),
            ('self', ['losses']),
            ('filters.losses_state_name', None),
            # ('self', ['debugging']),
            # ('filters.deltaG_fit_id', None),
            # ('filters.coverage_mse_fit_id', None),
        ]

    @param.depends('spin', watch=True)
    def _spin_updated(self):
        view = self.views['protein']
        view.ngl_view.spin = self.spin


class FileExportControl(ControlPanel):
    # todo check if docstring is true
    """
    This controller allows users to export and download datasets.

    All datasets can be exported as .txt tables.
    'Mappable' datasets (with r_number column) can be exported as .pml pymol script, which colors protein structures
    based on their 'color' column.

    """

    header = "File Export"
    table = param.Selector(label='Target dataset', doc='Name of the dataset to export')
    export_format = param.Selector(default='csv', objects=['csv', 'pprint'],
                                   doc="Format of the exported tables."
                                       "'csv' is machine-readable, 'pprint' is human-readable format")
    #todo add color param an dlink with protein viewer color

    def __init__(self, parent, **param):
        super(FileExportControl, self).__init__(parent, **param)

        objects = list(self.sources['dataframe'].tables.keys())
        self.param['table'].objects = objects
        self.table = objects[0]
        self.sources['dataframe'].param.watch(self._source_updated, 'updated')

    def make_dict(self):
        widgets = self.generate_widgets()
        widgets['export_tables'] = pn.widgets.FileDownload(
            label='Download table',
            callback=self.table_export_callback
        )
        widgets['export_pml'] = pn.widgets.FileDownload(label='Download pml scripts',
                                                        callback=self.pml_export_callback,
                                                        )

        return widgets

    @property
    def _layout(self):
        return [
            ('self', None)
        ]

    def _source_updated(self, *events):
        self.param['table'].objects = list(self.sources['dataframe'].tables.keys())
        self._table_updated()

    @param.depends('table', 'export_format', watch=True)
    def _table_updated(self):
        self.df = self.sources['dataframe'].get(self.table)

        ext = '.csv' if self.export_format == 'csv' else '.txt'
        self.widgets['export_tables'].filename = self.table + ext

        if self.table == 'colors':
            self.widgets['export_pml'].disabled = False
            self.widgets['export_pml'].filename = self.table + '_pml_scripts.zip'
        else:
            self.widgets['export_pml'].disabled = True

    @pn.depends('table')
    def pml_export_callback(self):

        if self.table:
            #todo check if table is valid for pml conversion

            bio = BytesIO()
            with zipfile.ZipFile(bio, 'w') as pml_zip:
                for col_name in self.df.columns:
                    name = col_name if isinstance(col_name, str) else '_'.join(col_name)
                    colors = self.df[col_name]
                    pml_script = series_to_pymol(colors)  # todo refactor pd_series_to_pymol?
                    pml_zip.writestr(name + '.pml', pml_script)

            bio.seek(0)
            return bio

    @pn.depends('table')  # param.depends?
    def table_export_callback(self):
        if self.table:
            io = dataframe_to_stringio(self.df, fmt=self.export_format)
            return io
        else:
            return None


class SingleMappingFileInputControl(MappingFileInputControl):
    """
    This controller allows users to upload *.txt files where quantities (protection factors, Gibbs free energy, etc) are
    mapped to a linear sequence.

    The column should be tab separated with on the last header line (starts with '#') the names of the columns. Columns
    should be tab-delimited.
    """

    def _action_add_dataset(self):
        super()._action_add_dataset()
        to_add_keys = set(self.parent.datasets.keys()) - set(self.parent.sources.keys())
        for key in to_add_keys:
            records = self.parent.datasets[key].to_records()
            data_source = DataSource(records, tags=['comparison', 'mapping'], x='r_number',
                                     renderer='circle', size=10)
            self.parent.publish_data(key, data_source)


class MatrixMappingFileInputControl(SingleMappingFileInputControl):
    datapoints = param.ListSelector(doc='Select datapoints to include in the matrix')

    def _action_add_dataset(self):
        super()._action_add_dataset()

        N = 20
        img = np.empty((N, N), dtype=np.uint32)
        view = img.view(dtype=np.uint8).reshape((N, N, 4))
        for i in range(N):
            for j in range(N):
                view[i, j, 0] = int(i / N * 255)
                view[i, j, 1] = 158
                view[i, j, 2] = int(j / N * 255)
                view[i, j, 3] = 255

        values = np.random.random(img.shape)

        img_ds_dict = {'img': [img], 'scores': [values]}
        data_source = DataSource(img_ds_dict, tags=['image'], name='scores_image', x=0, y=0)

        self.parent.publish_data('scores_image', data_source)

    def make_list(self):
        widget_list = super().make_list()
        datapoints_widget = widget_list.pop()
        widget_list.insert(3, datapoints_widget)
        return widget_list

    def _add_dataset(self):
        full_dict = self.protein.to_dict()
        data_dict = {k: v for k, v in full_dict.items() if k in self.datapoints}
        data_dict['r_number'] = self.protein.index
        protein = Protein(data_dict, index='r_number')
        self.parent.datasets[self.dataset_name] = protein

    @param.depends('input_file', watch=True)
    def _input_file_updated(self):
        super()._input_file_updated()
        if self.input_file:
            header_fields = self.protein.df.columns

            float_fields = [f for f in header_fields if f.replace('.', '', 1).isdigit()]
            self.param['datapoints'].objects = float_fields
            self.datapoints = float_fields

#        self.dataset_name = self.dataset_name or Path(self.widget_dict['input_file'].filename).stem


class MatrixImageControl(ControlPanel):
    """
    This controller takes an input loaded matrix and converts it to an (rgba) interpolated rendered image

    """


class FDPeptideFileInputControl(PeptideFileInputControl):
    # todo @tejas: Add test
    # This requires making a test function with the full_deuteration_app in apps.py
    def make_list(self):
        parameters = ['add_button', 'clear_button', 'drop_first', 'load_button', 'd_percentage',
                      'fd_state', 'fd_exposure', 'parse_button']
        first_widgets = list([self.widget_dict[par] for par in parameters])
        return self.file_selectors + first_widgets

    def _action_parse(self):
        """Apply controls to :class:`~pyhdx.models.PeptideMasterTable` and set :class:`~pyhdx.models.HDXMeasurement`"""
        pmt = self.parent.peptides

        data_states = pmt.data[pmt.data['state'] == self.fd_state]
        data_exposure = data_states[data_states['exposure'] == self.fd_exposure]

        scores = 100 * data_exposure['uptake'] / data_exposure['ex_residues']
        data_final = append_fields(data_exposure, 'scores', data=scores, usemask=False)

        # pmt.set_control((fd_state, fd_exposure))
        series = HDXMeasurement(data_final)

        self.parent.series = series

        self.parent.logger.info(f"Loaded FD control '{self.exp_state}' with {len(series.coverage)} peptides")
        self.parent.logger.info(f'Mean deuteration is {scores.mean()}%, std {scores.std()}%')


class PeptideFoldingFileInputControl(PeptideFileInputControl):
    # todo @tejas: Add test
    # This requires making a test function with the folding in apps.py

    be_mode = param.Selector(doc='Select method of normalization', label='Norm mode', objects=['Exp', 'Theory']
                             , precedence=-1)
    fd_state = param.Selector(doc='State used to normalize uptake', label='100% Control State')
    fd_exposure = param.Selector(doc='Exposure used to normalize uptake', label='100% Control Exposure')
    zero_state = param.Selector(doc='State used to zero uptake', label='0% Control State')
    zero_exposure = param.Selector(doc='Exposure used to zero uptake', label='0% Control Exposure')

    def make_dict(self):
        return self.generate_widgets()

    def make_list(self):
        parameters = ['add_button', 'clear_button', 'drop_first', 'ignore_prolines', 'load_button',
                      'fd_state', 'fd_exposure', 'zero_state', 'zero_exposure', 'exp_state',
                      'exp_exposures', 'parse_button']
        first_widgets = list([self.widget_dict[par] for par in parameters])
        return self.file_selectors + first_widgets

    def _action_load(self):
        super()._action_load()
        states = list(np.unique(self.parent.peptides.data['state']))
        self.param['zero_state'].objects = states
        self.zero_state = states[0]

    @param.depends('fd_state', 'fd_exposure', watch=True)
    def _update_experiment(self):
        #TODO THIS needs to be updated to also incorporate the zero (?)
        pm_dict = self.parent.peptides.return_by_name(self.fd_state, self.fd_exposure)
        states = list(np.unique([v.state for v in pm_dict.values()]))
        self.param['exp_state'].objects = states
        self.exp_state = states[0] if not self.exp_state else self.exp_state

    @param.depends('zero_state', watch=True)
    def _update_zero_exposure(self):
        b = self.parent.peptides.data['state'] == self.zero_state
        data = self.parent.peptides.data[b]
        exposures = list(np.unique(data['exposure']))
        self.param['zero_exposure'].objects = exposures
        if exposures:
            self.control_exposure = exposures[0]

    def _action_parse(self):
        """Apply controls to :class:`~pyhdx.models.PeptideMasterTable` and set :class:`~pyhdx.models.HDXMeasurement`"""
        control_0 = self.zero_state, self.zero_exposure
        self.parent.peptides.set_control((self.fd_state, self.fd_exposure), control_0=control_0)

        data_states = self.parent.peptides.data[self.parent.peptides.data['state'] == self.exp_state]
        data = data_states[np.isin(data_states['exposure'], self.exp_exposures)]

        series = HDXMeasurement(data)
        self.parent.series = series

        self._publish_scores()

        self.parent.logger.info(f'Loaded experiment state {self.exp_state} '
                                f'({len(series)} timepoints, {len(series.coverage)} peptides each)')


class DifferenceControl(ControlPanel):
    """
    This controller allows users to select two datasets from available datasets, choose a quantity to compare between,
    and choose the type of operation between quantities (Subtract/Divide).

    """
    header = 'Differences'

    dataset_1 = param.Selector(doc='First dataset to compare')
    dataset_2 = param.Selector(doc='Second dataset to compare')

    comparison_name = param.String()
    operation = param.Selector(default='Subtract', objects=['Subtract', 'Divide'],
                               doc='Select the operation to perform between the two datasets')

    comparison_quantity = param.Selector(doc="Select a quantity to compare (column from input txt file)")
    add_comparison = param.Action(lambda self: self._action_add_comparison(),
                                  doc='Click to add this comparison to available comparisons')
    comparison_list = param.ListSelector(doc='Lists available comparisons')
    remove_comparison = param.Action(lambda self: self._action_remove_comparison(),
                                     doc='Remove selected comparisons from the list')

    def __init__(self, parent, **params):
        super(DifferenceControl, self).__init__(parent, **params)
        self.parent.param.watch(self._datasets_updated, ['datasets'])

    def _datasets_updated(self, events):
        objects = list(self.parent.datasets.keys())

        self.param['dataset_1'].objects = objects
        if not self.dataset_1:
            self.dataset_1 = objects[0]
        self.param['dataset_2'].objects = objects
        if not self.dataset_2:# or self.dataset_2 == objects[0]:  # dataset2 default to second dataset? toggle user modify?
            self.dataset_2 = objects[0]

    @param.depends('dataset_1', 'dataset_2', watch=True)
    def _selection_updated(self):
        if self.datasets:
            unique_names = set.intersection(*[{name for name in protein.df.dtypes.index} for protein in self.datasets])
            objects = [name for name in unique_names if np.issubdtype(self.protein_1[name].dtype, np.number)]
            objects.sort()

            # todo check for scara dtype
            self.param['comparison_quantity'].objects = objects
            if self.comparison_quantity is None:
                self.comparison_quantity = objects[0]

    @property
    def protein_1(self):
        """:class:`~pyhdx.models.Protein`: Protein object of dataset 1"""
        try:
            return self.parent.datasets[self.dataset_1]
        except KeyError:
            return None

    @property
    def protein_2(self):
        """:class:`~pyhdx.models.Protein`: Protein object of dataset 2"""
        try:
            return self.parent.datasets[self.dataset_2]
        except KeyError:
            return None

    @property
    def datasets(self):
        """:obj:`tuple`: Tuple with `(protein_1, protein_2)"""
        datasets = (self.protein_1, self.protein_2)
        if None in datasets:
            return None
        else:
            return datasets

    def _action_add_comparison(self):
        if not self.comparison_name:
            self.parent.logger.info('The added comparison needs to have a name')
            return
        if self.datasets is None:
            return

        op = {'Subtract': operator.sub, 'Divide': operator.truediv}[self.operation]
        comparison = op(*[p[self.comparison_quantity] for p in self.datasets]).rename('comparison')
        value1 = self.protein_1[self.comparison_quantity].rename('value1')
        value2 = self.protein_2[self.comparison_quantity].rename('value2')
        df = pd.concat([comparison, value1, value2], axis=1)

        output = df.to_records()
        data_source = DataSource(output, tags=['comparison', 'mapping'], x='r_number', y='comparison',
                                 renderer='circle', size=10)
        self.parent.publish_data(self.comparison_name, data_source)  # Triggers parent.sources param
        self.comparison_name = ''

    def _action_remove_comparison(self):
        for comparison in self.comparison_list:
            self.parent.sources.pop(comparison)   #Popping from dicts does not trigger param
        self.parent.param.trigger('sources')

    @param.depends('parent.sources', watch=True)
    def _update_comparison_list(self):
        objects = [name for name, d in self.parent.sources.items() if 'comparison' in d.tags]
        self.param['comparison_list'].objects = objects


class SingleControl(ControlPanel):
    # todo @tejas: Add test

    """
    This controller allows users to select a dataset from available datasets, and choose a quantity to classify/visualize,
    and add this quantity to the available datasets.
    """

    #todo subclass with DifferenceControl
    #rename dataset_name
    header = 'Datasets'

    dataset = param.Selector(doc='Dataset')
    dataset_name = param.String(doc='Name of the dataset to add')
    quantity = param.Selector(doc="Select a quantity to plot (column from input txt file)")

    add_dataset = param.Action(lambda self: self._action_add_dataset(),
                               doc='Click to add this comparison to available comparisons')
    dataset_list = param.ListSelector(doc='Lists available comparisons')
    remove_dataset = param.Action(lambda self: self._action_remove_comparison(),
                                  doc='Remove selected datasets from available datasets')

    def __init__(self, parent, **params):
        super(SingleControl, self).__init__(parent, **params)
        self.parent.param.watch(self._datasets_updated, ['datasets'])

    def _datasets_updated(self, events):
        objects = list(self.parent.datasets.keys())

        self.param['dataset'].objects = objects
        if not self.dataset:
            self.dataset = objects[0]

    @param.depends('dataset', watch=True)
    def _selection_updated(self):
        if self.dataset:
            dataset = self.parent.datasets[self.dataset]
            names = dataset.dtype.names
            objects = [name for name in names if name != 'r_number']
            self.param['quantity'].objects = objects
            if self.quantity is None:
                self.quantity = objects[0]

    def _action_add_dataset(self):
        if not self.dataset_name:
            self.parent.logger.info('The added comparison needs to have a name')
            return
        if not self.dataset:
            return

        array = self.parent.datasets[self.dataset]
        data_source = DataSource(array, tags=['comparison', 'mapping'], x='r_number', y=self.quantity,
                                 renderer='circle', size=10)
        self.parent.publish_data(self.dataset_name, data_source)  # Triggers parent.sources param
        self.comparison_name = ''

    def _action_remove_comparison(self):
        for ds in self.dataset_list:
            self.parent.sources.pop(ds)   #Popping from dicts does not trigger param
        self.parent.param.trigger('sources')

    @param.depends('parent.sources', watch=True)
    def _update_dataset_list(self):
        objects = [name for name, d in self.parent.sources.items()]
        self.param['dataset_list'].objects = objects


class FDCoverageControl(CoverageControl):
    def make_list(self):
        lst = super(CoverageControl, self).make_list()
        return lst[:-1]


class FoldingFitting(InitialGuessControl):
    fitting_model = param.Selector(default='Dissociation', objects=['Dissociation'],
                                   doc='Choose method for determining initial guesses.')

    def make_list(self):
        self.widget_dict.update(pbar1=self.pbar1.view, pbar2=self.pbar2.view)
        parameters = ['fitting_model', 'lower_bound', 'upper_bound', 'do_fit1', 'pbar1']

        widget_list = list([self.widget_dict[par] for par in parameters])
        return widget_list


class FitResultControl(ControlPanel):
    # @tejas skip test, currently bugged, issue #182

    """
    This controller allows users to view to fit result and how it describes the uptake of every peptide.
    """

    header = 'Fit Results'

    peptide_index = param.Integer(0, bounds=(0, None),
                                 doc='Index of the peptide to display.')
    x_axis_type = param.Selector(default='Log', objects=['Linear', 'Log'],
                                 doc='Choose whether to plot the x axis as Logarithmic axis or Linear.')

    def __init__(self, parent, **param):
        super(FitResultControl, self).__init__(parent, **param)

        self.d_uptake = {}  ## Dictionary of arrays (N_p, N_t) with results of fit result model calls
        #todo why does still still exists should it not just be dataobjects??
        # --> because they need to be calcualted only once and then dataobjects are generated per index
        # can be improved probably (by putting all data in data source a priory?

        self.parent.param.watch(self._series_updated, ['datasets']) #todo refactor
        self.parent.param.watch(self._fit_results_updated, ['fit_results'])

    def _series_updated(self, *events):
        pass
        #
        # self.param['peptide_index'].bounds = (0, len(self.parent.series.coverage.data) - 1)
        # self.d_uptake['uptake_corrected'] = self.parent.series.uptake_corrected.T
        # self._update_sources()

    @property
    def fit_timepoints(self):
        time = np.logspace(-2, np.log10(self.parent.series.timepoints.max()), num=250)
        time = np.insert(time, 0, 0.)
        return time

    def _fit_results_updated(self, *events):
        accepted_fitresults = ['fr_pfact']
        #todo wrappertje which checks with a cached previous version of this particular param what the changes are even it a manual trigger
        for name, fit_result in self.parent.fit_results.items():
            if name in accepted_fitresults:
                D_upt = fit_result(self.fit_timepoints)
                self.d_uptake[name] = D_upt
            else:
                continue
        # push results to graph
            self._update_sources()

    @param.depends('peptide_index', watch=True)
    def _update_sources(self):
        for name, array in self.d_uptake.items():
            if name == 'uptake_corrected':  ## this is the raw data
                timepoints = self.parent.series.timepoints
                renderer = 'circle'
                color = '#000000'
            else:
                timepoints = self.fit_timepoints
                renderer = 'line'
                color = '#bd0d1f'  #todo css / default color cycle per Figure Panel?

            dic = {'time': timepoints, 'uptake': array[self.peptide_index, :]}
            data_source = DataSource(dic, x='time', y='uptake', tags=['uptake_curve'], renderer=renderer, color=color)
            self.parent.publish_data(name, data_source)


class ColoringControl(ClassificationControl):
    # WIP class, skip tests


    def make_dict(self):
        widgets_dict = super().make_dict()
        widgets_dict.pop('quantity')

        return widgets_dict

    @param.depends('values', 'colors', 'target', 'quantity', watch=True)
    def _get_colors(self):
        # todo this part is repeated
        if np.all(self.values == 0):
            return
        elif np.any(np.diff(self.values) > 0):  # Skip applying colors when not strictly monotonic descending
            return
        elif not self.target:
            return
        elif 'scores_image' not in self.parent.sources.keys():
            return

        tgt_source = self.parent.sources[self.target] # full array including nan entries
        r_number = tgt_source.source.data['r_number']
        assert np.all(np.diff(r_number) == 1)


        headers = [f for f in tgt_source.source.data.keys() if f.replace('.', '', 1).isdigit()]

        headers.sort(key=float)
        timepoints = np.array([float(f) for f in headers])
        N_interpolate = 500
        interp_timepoints = np.linspace(0, timepoints.max(), num=N_interpolate, endpoint=True)
        data_array = np.stack([tgt_source.source.data[k] for k in headers])

        array = np.stack([np.interp(interp_timepoints, timepoints, data) for data in data_array.T]).T


        colors_hex = self._calc_colors(array.flatten())  # colors are in hex format
        if colors_hex is None:  # this is the colors not between 0 and 1 bug / error
            return

        colors_hex[colors_hex == 'nan'] = '#8c8c8c'
        colors_rgba = np.array([hex_to_rgba(h) for h in colors_hex])

        shape = (N_interpolate, len(r_number))
        img = np.empty(shape, dtype=np.uint32)
        view = img.view(dtype=np.uint8).reshape(*shape, 4)
        view[:] = colors_rgba.reshape(*shape, 4)

        img_source = self.parent.sources['scores_image']
        img_source.render_kwargs['dw'] = r_number.max()
        img_source.render_kwargs['dh'] = timepoints.max()
        img_source.source.data.update(img=[img], scores=[array])


        #self.parent.sources[self.target].source.data['color'] = colors


class DifferenceFileExportControl(FileExportControl):
    """
    This controller allows users to export and download datasets.

    'Mappable' datasets (with r_number column) can be exported as .pml pymol script, which colors protein structures
    based on their 'color' column.

    """

    accepted_tags = ['mapping']
    #todo include comparison info (x vs y) in output

    def _sources_updated(self, *events):  #refactor _parent_sources_updated on classificationcontrol
        data_sources = [k for k, src in self.parent.sources.items() if src.resolve_tags(self.accepted_tags)]
        self.param['target'].objects = list(data_sources)

        # Set target if its not set already
        if not self.target and data_sources:
            self.target = data_sources[-1]

    @pn.depends('target', watch=True)
    def _update_filename(self):
        self.export_linear_download.filename = self.target + '_linear.txt'
        if 'r_number' in self.export_dict.keys():
            self.pml_script_download.filename = self.target + '_pymol.pml'


class OptionsControl(ControlPanel):
    """The controller is used for various settings."""

    header = 'Options'

    #todo this should be a component (mixin?) for apps who dont have these figures
    link_xrange = param.Boolean(True, doc='Link the X range of the coverage figure and other linear mapping figures.', constant=False)
    log_level = param.Selector(default='DEBUG', objects=['DEBUG', 'INFO', 'WARN', 'ERROR', 'FATAL', 'OFF', 'TRACE'],
                               doc='Set the logging level.')

    def __init__(self, parent, **param):
        super(OptionsControl, self).__init__(parent, **param)

    @property
    def enabled(self):
        return self.master_figure is not None and self.client_figures is not None

    @param.depends('link_xrange', watch=True)
    def _update_link(self):
        if self.enabled:
            if self.link_xrange:
                self._link()
            else:
                self._unlink()

    @property
    def client_figures(self):
        client_names = ['RateFigure', 'PFactFigure']
        return [self.parent.figure_panels[name].figure for name in client_names]

    @property
    def master_figure(self):
        return self.parent.figure_panels['CoverageFigure'].figure

    @property
    def figures(self):
        return [self.master_figure] + self.client_figures

    def _unlink(self):
        for fig in self.figures:
            fig.x_range.js_property_callbacks.pop('change:start')
            fig.x_range.js_property_callbacks.pop('change:end')

    def _link(self):
        for client in self.client_figures:
            self.master_figure.x_range.js_link('start',  client.x_range, 'start')
            self.master_figure.x_range.js_link('end', client.x_range, 'end')

            client.x_range.js_link('start', self.master_figure.x_range, 'start')
            client.x_range.js_link('end', self.master_figure.x_range, 'end')


class DeveloperControl(ControlPanel):
    """Controller with debugging options"""

    header = 'Developer Options'
    test_logging = param.Action(lambda self: self._action_test_logging())
    breakpoint_btn = param.Action(lambda self: self._action_break())
    test_btn = param.Action(lambda self: self._action_test())
    trigger_btn = param.Action(lambda self: self._action_trigger())
    print_btn = param.Action(lambda self: self._action_print())
    runtime_warning = param.Action(lambda self: self._action_runtime())

    def __init__(self, parent, **params):
        super(DeveloperControl, self).__init__(parent, **params)

    def _action_test_logging(self):
        print(self.parent.logger)
        self.parent.logger.debug('TEST DEBUG MESSAGE')
        #logging.info('THis is some info')
        for i in range(20):
            self.parent.logger.info('dit is een test123')

    def _action_print(self):

        hdx_set = self.parent.hdx_set
        print(hdx_set.names)
        guess = self.parent.control_panels['FitControl']
        rates_df = self.sources['dataframe'].get('rates', fit_ID=guess.initial_guess)
        print(guess.initial_guess)
        print(rates_df)

        rates_guess = [rates_df[state]['rate'] for state in hdx_set.names]
        gibbs_guess = hdx_set.guess_deltaG(rates_guess)

    def _action_break(self):
        main_ctrl = self.parent
        control_panels = main_ctrl.control_panels
        views = main_ctrl.views
        sources = main_ctrl.sources

        mse_view = views['coverage_mse']
        data = mse_view.get_data()
        print('mse')
        print(data)

        coverage_view = views['coverage']
        data = coverage_view.get_data()
        print('coverage')
        print(data)


        print('Time for a break')

    def _action_test(self):
        src_file = r'C:\Users\jhsmi\pp\PyHDX\tests\test_data\ecSecB_torch_fit.txt'
        array = txt_to_np(src_file)
        data_dict = {name: array[name] for name in array.dtype.names}

        data_dict['color'] = np.full_like(array, fill_value=DEFAULT_COLORS['pfact'], dtype='<U7')
        data_source = DataSource(data_dict, x='r_number', tags=['mapping', 'pfact', 'deltaG'],
                                 renderer='circle', size=10, name='global_fit')

        self.parent.publish_data('global_fit', data_source)

    def _action_trigger(self):
        deltaG_figure = self.parent.figure_panels['DeltaGFigure']
        deltaG_figure.bk_pane.param.trigger('object')

    def _action_runtime(self):
        result = np.mean([])