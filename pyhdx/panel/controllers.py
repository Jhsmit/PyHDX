from pyhdx.models import PeptideMasterTable, KineticsSeries, Protein
from pyhdx.panel.widgets import NumericInput
from pyhdx.panel.data_sources import DataSource
from pyhdx.panel.base import ControlPanel, DEFAULT_COLORS, DEFAULT_CLASS_COLORS
from pyhdx.fitting import KineticsFitting
from pyhdx.fileIO import read_dynamx, txt_to_np, fmt_export
from pyhdx.support import autowrap, colors_to_pymol, rgb_to_hex
from pyhdx import VERSION_STRING
from scipy import constants
import param
import panel as pn
import numpy as np
from numpy.lib.recfunctions import append_fields
from pathlib import Path
from skimage.filters import threshold_multiotsu
from numpy.lib.recfunctions import stack_arrays
from io import StringIO
from tornado.ioloop import IOLoop
from functools import partial
from bokeh.models import ColumnDataSource, LinearColorMapper, ColorBar
from bokeh.plotting import figure
from collections import namedtuple
import operator
import matplotlib as mpl
import matplotlib.pyplot as plt
import itertools
import pandas as pd

from .widgets import ColoredStaticText, ASyncProgressBar

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
        self.dataset_name = self.dataset_name or Path(self._widget_dict['input_file'].filename).stem

    #todo refactor dataset to protein_something
    def _action_add_dataset(self):
        if self.dataset_name in self.parent.datasets.keys():
            self.parent.logger.info(f'Dataset {self.dataset_name} already added')
        elif not self.dataset_name:
            self.parent.logger.info('The added comparison needs to have a name')
        elif not self.input_file:
            self.parent.logger.info('Empty or no file selected')
        else:
            try:
                sio = StringIO(self.input_file.decode())
                array = txt_to_np(sio)
                assert 'r_number' in array.dtype.names, "Input file needs to have an 'r_number' column"
                array['r_number'] += self.offset
                protein = Protein(array, index='r_number')
                self.parent.datasets[self.dataset_name] = protein
                self.parent.param.trigger('datasets')
            except UnicodeDecodeError:
                self.parent.logger.info('Invalid file type, supplied file is not a text file')

        self._widget_dict['input_file'].filename = ''
        self._widget_dict['input_file'].value = b''

        self.dataset_name = ''

    def _action_remove_dataset(self):
        for dataset_name in self.datasets_list:
            self.parent.datasets.pop(dataset_name)
        self.parent.param.trigger('datasets')

    def _datasets_updated(self, events):
        self.param['datasets_list'].objects = list(self.parent.datasets.keys())


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


class PeptideFileInputControl(ControlPanel):
    """
    This controller allows users to input .csv file (Currently only DynamX format) of 'state' peptide uptake data.
    Users can then choose how to correct for back-exchange and which 'state' and exposure times should be used for
    analysis.

    """
    header = 'Peptide Input'

    add_button = param.Action(lambda self: self._action_add(), doc='Add File', label='Add File')
    clear_button = param.Action(lambda self: self._action_clear(), doc='Clear files', label='Clear Files')
    drop_first = param.Integer(1, bounds=(0, None), doc='Select the number of N-terminal residues to ignore.')
    ignore_prolines = param.Boolean(True, constant=True, doc='Prolines are ignored as they do not exchange D.')
    d_percentage = param.Number(95., bounds=(0, 100), doc='Percentage of deuterium in the labelling buffer',
                                label='Deuterium percentage')
    load_button = param.Action(lambda self: self._action_load(), doc='Load the selected files', label='Load Files')

    be_mode = param.Selector(doc='Select method of back exchange correction', label='Norm mode', objects=['Exp', 'Theory'])
    fd_state = param.Selector(doc='State used to normalize uptake', label='FD State')
    fd_exposure = param.Selector(doc='Exposure used to normalize uptake', label='FD Exposure')
    be_percent = param.Number(28., bounds=(0, 100), doc='Global percentage of back-exchange',
                              label='Back exchange percentage')

    exp_state = param.Selector(doc='State for selected experiment', label='Experiment State')
    exp_exposures = param.ListSelector(default=[], objects=[''], label='Experiment Exposures'
                                       , doc='Selected exposure time to use')

    c_term = param.Integer(0, bounds=(0, None),
                           doc='Index of the c terminal residue in the protein. Used for generating pymol export script')
    parse_button = param.Action(lambda self: self._action_parse(), label='Parse',
                                doc='Parse selected peptides for further analysis and apply back-exchange correction')

    def __init__(self, parent, **params):
        self.file_selectors = [pn.widgets.FileInput(accept='.csv')]
        super(PeptideFileInputControl, self).__init__(parent, **params)

    def make_dict(self):
        return self.generate_widgets(be_mode=pn.widgets.RadioButtonGroup, be_percent=pn.widgets.FloatInput,
                                     d_percentage=pn.widgets.FloatInput)

    def make_list(self):
        parameters = ['add_button', 'clear_button', 'drop_first', 'ignore_prolines', 'd_percentage', 'load_button',
                      'be_mode', 'fd_state', 'fd_exposure', 'exp_state',
                      'exp_exposures', 'c_term', 'parse_button']
        first_widgets = list([self._widget_dict[par] for par in parameters])
        return self.file_selectors + first_widgets

    def _action_add(self):
        """Add another FileInput widget/"""
        widget = pn.widgets.FileInput(accept='.csv')
        i = len(self.file_selectors)  # position to insert the new file selector into the widget box
        self.file_selectors.append(widget)
        self._box.insert(i, widget)

    def _action_clear(self):
        """Clear all file selectors and set number of file selectors to one."""
        self.parent.logger.debug('Cleared file selectors')

        while self.file_selectors:
            fs = self.file_selectors.pop()
            idx = list(self._box).index(fs)
            self._box.pop(idx)
        self._action_add()

    def _action_load(self):
        """Load files from FileInput widgets """
        data_list = []
        for file_selector in self.file_selectors:
            if file_selector.value is not None:
                s_io = StringIO(file_selector.value.decode('UTF-8'))
                data = read_dynamx(s_io)
                data_list.append(data)

        combined = stack_arrays(data_list, asrecarray=True, usemask=False, autoconvert=True)

        self.parent.peptides = PeptideMasterTable(combined, d_percentage=self.d_percentage,
                                                  drop_first=self.drop_first, ignore_prolines=self.ignore_prolines)

        states = list(np.unique(self.parent.peptides.data['state']))
        self.param['fd_state'].objects = states
        self.fd_state = states[0]
        #self.param['zero_state'].objects = ['None'] + states
        #self.zero_state = 'None'

        self.parent.logger.info(
            f'Loaded {len(data_list)} file{"s" if len(data_list) > 1 else ""} with a total '
            f'of {len(self.parent.peptides)} peptides')

    def _action_parse(self):
        """Apply controls to :class:`~pyhdx.models.PeptideMasterTable` and set :class:`~pyhdx.models.KineticsSeries`"""
        if self.be_mode == 'Exp':
            control_0 = None # = (self.zero_state, self.zero_exposure) if self.zero_state != 'None' else None
            self.parent.peptides.set_control((self.fd_state, self.fd_exposure), control_0=control_0)
        elif self.be_mode == 'Theory':
            self.parent.peptides.set_backexchange(self.be_percent)

        data_states = self.parent.peptides.data[self.parent.peptides.data['state'] == self.exp_state]
        data = data_states[np.isin(data_states['exposure'], self.exp_exposures)]

        series = KineticsSeries(data, c_term=self.c_term)
        self.parent.series = series

        self.parent.logger.info(f'Loaded experiment state {self.exp_state} '
                                f'({len(series)} timepoints, {len(series.cov)} peptides each)')

    @param.depends('be_mode', watch=True)
    def _update_be_mode(self):
        if self.be_mode == 'Exp':
            self.box_pop('be_percent')
            self.box_insert_after('be_mode', 'fd_state')
            self.box_insert_after('fd_state', 'fd_exposure')

        elif self.be_mode == 'Theory':
            self.box_pop('fd_state')
            self.box_pop('fd_exposure')
            self.box_insert_after('be_mode', 'be_percent')

            try:
                states = np.unique(self.parent.peptides.data['state'])
                self.param['exp_state'].objects = states
                self.exp_state = states[0] if not self.exp_state else self.exp_state
            except (TypeError, AttributeError):
                pass

    @param.depends('fd_state', watch=True)
    def _update_norm_exposure(self):
        b = self.parent.peptides.data['state'] == self.fd_state
        data = self.parent.peptides.data[b]
        exposures = list(np.unique(data['exposure']))
        self.param['fd_exposure'].objects = exposures
        if exposures:
            self.fd_exposure = exposures[0]

    #todo refactor norm to FD
    @param.depends('fd_state', 'fd_exposure', watch=True)
    def _update_experiment(self):
        #TODO THIS needs to be updated to also incorporate the zero (?)
        pm_dict = self.parent.peptides.return_by_name(self.fd_state, self.fd_exposure)
        states = list(np.unique([v.state for v in pm_dict.values()]))
        self.param['exp_state'].objects = states
        self.exp_state = states[0] if not self.exp_state else self.exp_state

    @param.depends('exp_state', watch=True)
    def _update_experiment_exposure(self):
        b = self.parent.peptides.data['state'] == self.exp_state
        exposures = list(np.unique(self.parent.peptides.data['exposure'][b]))
        exposures.sort()
        self.param['exp_exposures'].objects = exposures
        self.exp_exposures = exposures

        self.c_term = int(np.max(self.parent.peptides.data['end'][b]))


class FDPeptideFileInputControl(PeptideFileInputControl):

    def make_list(self):
        parameters = ['add_button', 'clear_button', 'drop_first', 'load_button', 'd_percentage',
                      'fd_state', 'fd_exposure', 'parse_button']
        first_widgets = list([self._widget_dict[par] for par in parameters])
        return self.file_selectors + first_widgets

    def _action_parse(self):
        """Apply controls to :class:`~pyhdx.models.PeptideMasterTable` and set :class:`~pyhdx.models.KineticsSeries`"""
        pmt = self.parent.peptides

        data_states = pmt.data[pmt.data['state'] == self.fd_state]
        data_exposure = data_states[data_states['exposure'] == self.fd_exposure]

        scores = 100 * data_exposure['uptake'] / data_exposure['ex_residues']
        data_final = append_fields(data_exposure, 'scores', data=scores, usemask=False)

        # pmt.set_control((fd_state, fd_exposure))
        series = KineticsSeries(data_final)

        self.parent.series = series

        self.parent.logger.info(f"Loaded FD control '{self.exp_state}' with {len(series.cov)} peptides")
        self.parent.logger.info(f'Mean deuteration is {scores.mean()}%, std {scores.std()}%')


class PeptideFoldingFileInputControl(PeptideFileInputControl):
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
        first_widgets = list([self._widget_dict[par] for par in parameters])
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
        """Apply controls to :class:`~pyhdx.models.PeptideMasterTable` and set :class:`~pyhdx.models.KineticsSeries`"""
        control_0 = self.zero_state, self.zero_exposure
        self.parent.peptides.set_control((self.fd_state, self.fd_exposure), control_0=control_0)

        data_states = self.parent.peptides.data[self.parent.peptides.data['state'] == self.exp_state]
        data = data_states[np.isin(data_states['exposure'], self.exp_exposures)]

        series = KineticsSeries(data)
        self.parent.series = series

        exposure_data_dict = {str(exposure): dpt for dpt, exposure in zip(series.scores_stack, series.timepoints)}
        data_dict = dict(r_number=series.cov.r_number, **exposure_data_dict)

        data_source = DataSource(data_dict, x='r_number', y=list(exposure_data_dict.keys()), tags=['mapping', 'scores'],
                                 renderer='circle', size=10, name='scores')
        self.parent.publish_data('scores', data_source)

        self.parent.logger.info(f'Loaded experiment state {self.exp_state} '
                                f'({len(series)} timepoints, {len(series.cov)} peptides each)')


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


class CoverageControl(ControlPanel):
    """
    This controller allows users to control the peptide coverage figure, by choosing how many peptides to plot vertically,
    which color map to use, and which exposure time to show.
    """
    header = 'Coverage'

    wrap = param.Integer(25, bounds=(0, None), doc='Number of peptides vertically before moving to the next row.') # todo auto?
    color_map = param.Selector(objects=['jet', 'inferno', 'viridis', 'cividis', 'plasma', 'cubehelix'], default='jet',
                               doc='Color map for coloring peptides by their deuteration percentage.')
    index = param.Integer(0, bounds=(0, 10), doc='Current index of coverage plot in time.')

    def __init__(self, parent, **params):
        self.exposure_str = ColoredStaticText(name='Exposure', value='0')  # todo update to some param?

        # We need a reference to color mapper to update it when the cmap changes
        self.color_mapper = LinearColorMapper(palette=self.palette, low=0, high=100)
        self.color_bar = self.get_color_bar()

        super(CoverageControl, self).__init__(parent, **params)
        self.parent.param.watch(self._series_updated, ['series'])

    def make_list(self):
        lst = super(CoverageControl, self).make_list()
        return lst + [self.exposure_str]#\, self.color_bar]

    def make_dict(self):
        return self.generate_widgets(index=pn.widgets.IntSlider)

    @param.depends('color_map', watch=True)
    def _update_cbar(self):
        cmap = mpl.cm.get_cmap(self.color_map)
        pal = tuple(mpl.colors.to_hex(cmap(value)) for value in np.linspace(0, 1, 1024, endpoint=True))
        self.color_mapper.palette = pal

    def get_color_bar(self):
        """pn.pane.Bokeh: bokeh pane with empty figure and only a color bar
        Currently is buggy when added in ``make_list``
        """
        # f5f5f5
        # from default value in panel css
        #https://github.com/holoviz/panel/blob/67bf192ea4138825ab9682c8f38bfe2d696a4e9b/panel/_styles/widgets.css
        color_bar = ColorBar(color_mapper=self.color_mapper, location=(0, 0), orientation='horizontal',
                             background_fill_color='#f5f5f5')
                             #title='D uptake percentage',  title_text_align='center')
        fig = figure(toolbar_location=None, title='D uptake percentage')
        fig.title.align = 'center'
        fig.background_fill_color = '#f5f5f5'
        fig.border_fill_color = '#f5f5f5'
        fig.outline_line_color = None
        fig.add_layout(color_bar, 'above')

        return pn.pane.Bokeh(fig, height=100, width=350)

    @property
    def palette(self):
        """`obj`:tuple: Tuple of hex colors"""
        cmap = mpl.cm.get_cmap(self.color_map)
        pal = tuple(mpl.colors.to_hex(cmap(value)) for value in np.linspace(0, 1, 1024, endpoint=True))
        return pal

    @property
    def peptide_measurement(self):
        if self.parent.series is not None:
            return self.parent.series[self.index]
        else:
            return None

    @property
    def coverage(self):
        """Coverage object describing the peptide layout"""
        return self.parent.series.cov

    @property
    def color(self):
        """~class:`np.ndarray`: array of color for each peptide based on their uptake score"""
        cmap = mpl.cm.get_cmap(self.color_map)
        c_rgba = cmap(self.peptide_measurement.data['scores'] / 100)
        c = [mpl.colors.to_hex(color) for color in c_rgba]

        return np.array(c)

    def _series_updated(self, event):
        """Triggered when there is a new :class:`~pyhdx.models.KineticsSeries loaded on the main controller""" #todo refactor
        # series must be uniform
        self.wrap = autowrap(self.coverage)
        self.param['index'].bounds = (0, len(event.new) - 1)

        # set index to zero
        self.index = 0

        width = self.coverage.data['end'] - self.coverage.data['start'] # Bars are inclusive, inclusive
        x = self.coverage.data['start'] - 0.5 + (width / 2)
        y = list(itertools.islice(itertools.cycle(range(self.wrap, 0, -1)), len(self.coverage)))
        index = [str(i) for i in range(len(self.coverage))]

        #plot_dict = dict(x=x, y=y, width=width, color=self.color, index=index)
        prop_dict = {name: self.peptide_measurement.data[name] for name in self.peptide_measurement.data.dtype.names}
        prop_dict.update(color=self.color, index=index, width=width, x=x, y=y)  # if the names are x and y no need to specify through render_kwargs
        source = DataSource(prop_dict, tags=['coverage'], name=f'coverage_{self.parent.series.state}')

        self.parent.publish_data('coverage', source)

    @param.depends('wrap', watch=True)
    def _update_wrap(self):
        try:
            y = list(itertools.islice(itertools.cycle(range(self.wrap, 0, -1)), len(self.coverage)))
            self.parent.sources['coverage'].source.data.update(y=y)
        except (KeyError, AttributeError):
            pass

    @param.depends('index', 'color_map', watch=True)
    def _update_colors(self):
        try:
            self.exposure_str.value = str(self.peptide_measurement.exposure)  # todo this should be an js_link?
            tooltip_fields = {field: self.peptide_measurement.data[field] for field in ['scores', 'uptake', 'uptake_corrected']}
            self.parent.sources['coverage'].source.data.update(color=self.color, **tooltip_fields)

        except (KeyError, AttributeError):
            pass


class FDCoverageControl(CoverageControl):
    def make_list(self):
        lst = super(CoverageControl, self).make_list()
        return lst[:-1]


class InitialGuessControl(ControlPanel):
    """
    This controller allows users to derive initial guesses for D-exchange rate from peptide uptake data.
    """

    #todo remove lambda symbol although its really really funny
    header = 'Initial Guesses'
    fitting_model = param.Selector(default='Half-life (λ)', objects=['Half-life (λ)', 'Association'],
                                   doc='Choose method for determining initial guesses.')

    lower_bound = param.Number(0., doc='Lower bound for association model fitting')
    upper_bound = param.Number(0., doc='Upper bound for association model fitting')
    do_fit1 = param.Action(lambda self: self._action_fit(), label='Do fitting', doc='Start initial guess fitting',
                           constant=True)

    def __init__(self, parent, **params):
        self.pbar1 = ASyncProgressBar()
        self.pbar2 = ASyncProgressBar()
        super(InitialGuessControl, self).__init__(parent, **params)
        self.parent.param.watch(self._parent_series_updated, ['series'])

    def make_dict(self):
        return self.generate_widgets(lower_bound=pn.widgets.LiteralInput, upper_bound=pn.widgets.LiteralInput)

    def make_list(self):
        self._widget_dict.update(pbar1=self.pbar1.view, pbar2=self.pbar2.view)
        parameters = ['fitting_model', 'do_fit1', 'pbar1']

        widget_list = list([self._widget_dict[par] for par in parameters])
        return widget_list

    @param.depends('fitting_model', watch=True)
    def _fitting_model_updated(self):
        if self.fitting_model == 'Half-life (λ)':
            self.box_pop('lower_bound')
            self.box_pop('upper_bound')
        elif self.fitting_model in ['Association', 'Dissociation']:
            self.box_insert_after('fitting_model', 'upper_bound')
            self.box_insert_after('fitting_model', 'lower_bound')

    def _parent_series_updated(self, events):
        if self.parent.series is not None:
            self.param['do_fit1'].constant = False

            kf = KineticsFitting(self.parent.series)
            self.lower_bound, self.upper_bound = kf.bounds

    async def _fit1_async(self):
        """Do fitting asynchronously on (remote) cluster"""
        kf = KineticsFitting(self.parent.series, cluster=self.parent.cluster, bounds=(self.lower_bound, self.upper_bound))
        fit_result = await kf.weighted_avg_fit_async(model_type=self.fitting_model.lower(), pbar=self.pbar1)
        self.parent.fit_results['fit1'] = fit_result

        output = fit_result.output.to_records('r_number')  # todo remove in between numpy step
        #todo duplicate code in fit - > method on parent?
        dic = {name: output[name] for name in output.dtype.names}
        dic['y'] = output['rate']  # entry y is by default used for plotting and thresholding
        dic['color'] = np.full_like(output, fill_value=DEFAULT_COLORS['fit1'], dtype='<U7')

        data_source = DataSource(dic, x='r_number', y='rate', tags=['mapping', 'rate'],
                                 renderer='circle', size=10)

        #trigger plot update
        callback = partial(self.parent.publish_data, 'fit1', data_source)
        self.parent.doc.add_next_tick_callback(callback)

        with pn.io.unlocked():
             #self.parent.publish_data('fit1', data_source)
             self.parent.param.trigger('fit_results')  #informs other fittings that initial guesses are now available
             self.pbar1.reset()
             self.param['do_fit1'].constant = False

    def _fit1(self):
        kf = KineticsFitting(self.parent.series, bounds=(self.lower_bound, self.upper_bound))
        fit_result = kf.weighted_avg_fit(model_type=self.fitting_model.lower(), pbar=self.pbar1)
        self.parent.fit_results['fit1'] = fit_result
        output = fit_result.output.to_records('r_number')  # todo remove in between numpy step
        dic = {name: output[name] for name in output.dtype.names}
        dic['y'] = output['rate']  # entry y is by default used for plotting and thresholding
        dic['color'] = np.full_like(output, fill_value=DEFAULT_COLORS['fit1'], dtype='<U7')

        data_source = DataSource(dic, x='r_number', y='rate', tags=['mapping', 'rate'],
                                 renderer='circle', size=10, name='fit1')

        self.parent.publish_data('fit1', data_source)  #todo refactor to force require setting name on DataSource Ojbect
        self.parent.param.trigger('fit_results')  # Informs TF fitting that now fit1 is available as initial guesses

        self.param['do_fit1'].constant = False
        self.pbar1.reset()

    def _action_fit(self):
        if self.parent.series is None:
            self.parent.logger.debug('No dataset loaded')
            return

        self.parent.logger.debug('Start initial guess fit')
        #todo context manager?
        self.param['do_fit1'].constant = True

        if self.fitting_model == 'Half-life (λ)':
            kf = KineticsFitting(self.parent.series)
            output = kf.weighted_avg_t50()
            fit_result = HalfLifeFitResult(output=output)
            array = output.to_records()
            dic = {name: array[name] for name in array.dtype.names}

            #dic['y'] = output['rate']  # entry y is by default used for plotting and thresholding
            #todo colors dont work (because DataSource init)
            dic['color'] = np.full_like(array, fill_value=DEFAULT_COLORS['half-life'], dtype='<U7')

            data_source = DataSource(dic, x='r_number', y='rate', tags=['mapping', 'rate'],
                                     renderer='circle', size=10, name='half-life')

            self.parent.publish_data('half-life', data_source)
            self.parent.fit_results['half-life'] = fit_result

            self.parent.param.trigger('fit_results')  # Informs TF fitting that now fit1 is available as initial guesses

            self.param['do_fit1'].constant = False
        else:

            if self.parent.cluster:
                self.parent._doc = pn.state.curdoc
                loop = IOLoop.current()
                loop.add_callback(self._fit1_async)
            else:
                self._fit1()


class FoldingFitting(InitialGuessControl):
    fitting_model = param.Selector(default='Dissociation', objects=['Dissociation'],
                                   doc='Choose method for determining initial guesses.')

    def make_list(self):
        self._widget_dict.update(pbar1=self.pbar1.view, pbar2=self.pbar2.view)
        parameters = ['fitting_model', 'lower_bound', 'upper_bound', 'do_fit1', 'pbar1']

        widget_list = list([self._widget_dict[par] for par in parameters])
        return widget_list


class FitControl(ControlPanel):
    """
    This controller allows users to execute TensorFlow fitting of the global data set.

    Currently, repeated fitting overrides the old result.
    """

    header = 'Fitting'
    initial_guess = param.Selector(doc='Name of dataset to use for initial guesses.')
    temperature = param.Number(293.15, doc='Deuterium labelling temperature in Kelvin')
    pH = param.Number(8., doc='Deuterium labelling pH', label='pH')

    stop_loss = param.Number(0.05, bounds=(0, None),
                             doc='Threshold loss difference below which to stop fitting.')
    stop_patience = param.Integer(50, bounds=(1, None),
                                  doc='Number of epochs where stop loss should be satisfied before stopping.')
    learning_rate = param.Number(10, bounds=(0, None),
                                 doc='Learning rate parameter for optimization.')
    momentum = param.Number(0.5, bounds=(0, None),
                            doc='Stochastic Gradient Descent momentum')
    nesterov = param.Boolean(True, doc='Use Nesterov type of momentum for SGD')
    epochs = param.Number(100000, bounds=(1, None),
                          doc='Maximum number of epochs (iterations.')
    regularizer = param.Number(2, bounds=(0, None), doc='Value for the regularizer.')
    do_fit = param.Action(lambda self: self._do_fitting(), constant=True, label='Do Fitting',
                          doc='Start global fitting')

    def __init__(self, parent, **params):
        self.pbar1 = ASyncProgressBar()
        super(FitControl, self).__init__(parent, **params)
        self.parent.param.watch(self._parent_fit_results_updated, ['fit_results'])

    def _parent_fit_results_updated(self, *events):
        possible_initial_guesses = ['half-life', 'fit1']
        objects = [name for name in possible_initial_guesses if name in self.parent.fit_results.keys()]
        if objects:
            self.param['do_fit'].constant = False

        self.param['initial_guess'].objects = objects
        if not self.initial_guess and objects:
            self.initial_guess = objects[0]

    def _do_fitting(self):
        self.param['do_fit'].constant = True
        self.parent.logger.debug('Start PyTorch fit')

        kf = KineticsFitting(self.parent.series, temperature=self.temperature, pH=self.pH)
        initial_result = self.parent.fit_results[self.initial_guess].output   #todo initial guesses could be derived from the CDS rather than fit results object
        result = kf.global_fit_torch(initial_result, reg=self.regularizer, learning_rate=self.learning_rate,
                                     momentum=self.momentum, nesterov=self.nesterov, epochs=self.epochs,
                                     patience=self.stop_patience, stop_loss=self.stop_loss)

        output = result.output
        output_name = 'global_fit'
        output.df['color'] = np.full_like(output, fill_value=DEFAULT_COLORS['pfact'], dtype='<U7') #todo change how default colors are determined
        data_source = DataSource(output, x='r_number', tags=['mapping', 'pfact', 'deltaG'], name=output_name,
                                 renderer='circle', size=10)

        self.parent.fit_results['fr_' + output_name] = result
        self.parent.publish_data(output_name, data_source)

        self.param['do_fit'].constant = False
        self.parent.param.trigger('fit_results')

        self.parent.logger.debug('Finished PyTorch fit')
        loss = result.metadata['mse_loss']
        self.parent.logger.info(f"Finished fitting in {len(loss)} epochs, final mean squared residuals is {loss[-1]}")


class FitResultControl(ControlPanel):
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

        self.parent.param.watch(self._series_updated, ['series'])
        self.parent.param.watch(self._fit_results_updated, ['fit_results'])

    def _series_updated(self, *events):
        self.param['peptide_index'].bounds = (0, len(self.parent.series.cov.data) - 1)
        self.d_uptake['uptake_corrected'] = self.parent.series.uptake_corrected.T
        self._update_sources()

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


class ClassificationControl(ControlPanel):
    """
    This controller allows users classify 'mapping' datasets and assign them colors.

    Coloring can be either in discrete categories or as a continuous custom color map.
    """

    header = 'Classification'
    # format ['tag1', ('tag2a', 'tag2b') ] = tag1 OR (tag2a AND tag2b)
    accepted_tags = ['mapping']

    # todo unify name for target field (target_data set)
    # When coupling param with the same name together there should be an option to exclude this behaviour
    target = param.Selector(label='Target')
    quantity = param.Selector(label='Quantity')

    mode = param.Selector(default='Discrete', objects=['Discrete', 'Continuous'],
                          doc='Choose color mode (interpolation between selected colors).')#, 'ColorMap'])
    num_colors = param.Integer(3, bounds=(1, 10),
                              doc='Number of classification colors.')
    otsu_thd = param.Action(lambda self: self._action_otsu(), label='Otsu',
                            doc="Automatically perform thresholding based on Otsu's method.")
    linear_thd = param.Action(lambda self: self._action_linear(), label='Linear',
                              doc='Automatically perform thresholding by creating equally spaced sections.')
    log_space = param.Boolean(True,
                              doc='Boolean to set whether to apply colors in log space or not.')

    show_thds = param.Boolean(True, label='Show Thresholds', doc='Toggle to show/hide threshold lines.')
    values = param.List(precedence=-1)
    colors = param.List(precedence=-1)

    def __init__(self, parent, **param):
        super(ClassificationControl, self).__init__(parent, **param)

        self.values_widgets = []
        self.colors_widgets = []
        self._update_num_colors()
        self._update_num_values()

        self.param.trigger('values')
        self.param.trigger('colors')
        self.parent.param.watch(self._parent_sources_updated, ['sources'])

    def make_dict(self):
        return self.generate_widgets(num_colors=pn.widgets.IntInput, mode=pn.widgets.RadioButtonGroup)

    def _parent_sources_updated(self, *events):
        data_sources = [k for k, src in self.parent.sources.items() if src.resolve_tags(self.accepted_tags)]
        self.param['target'].objects = list(data_sources)

        # Set target if its not set already
        if not self.target and data_sources:
            self.target = data_sources[-1]

        if self.values:
            self._get_colors()

    @param.depends('target', watch=True)
    def _target_updated(self):
        data_source = self.parent.sources[self.target]
        self.param['quantity'].objects = [f for f in data_source.scalar_fields if not f.startswith('_')]
        default_priority = ['deltaG', 'comparison']  # Select these fields by default if they are present
        if not self.quantity and data_source.scalar_fields:
            for field in default_priority:
                if field in data_source.scalar_fields:
                    self.quantity = field
                    break
                self.quantity = data_source.scalar_fields[0]

    @property
    def target_array(self):
        """returns the array to calculate colors from, NaN entries are removed"""

        try:
            y_vals = self.parent.sources[self.target][self.quantity]
            return y_vals[~np.isnan(y_vals)]
        except KeyError:
            return None

    def _action_otsu(self):
        if self.num_colors > 1 and self.target_array is not None:
            func = np.log if self.log_space else lambda x: x  # this can have NaN when in log space
            thds = threshold_multiotsu(func(self.target_array), classes=self.num_colors)
            for thd, widget in zip(thds[::-1], self.values_widgets):  # Values from high to low
                widget.value = np.exp(thd) if self.log_space else thd
        self._get_colors()

    def _action_linear(self):
        i = 1 if self.mode == 'Discrete' else 0
        if self.log_space:
            thds = np.logspace(np.log(np.min(self.target_array)), np.log(np.max(self.target_array)),
                               num=self.num_colors + i, endpoint=True, base=np.e)
        else:
            thds = np.linspace(np.min(self.target_array), np.max(self.target_array), num=self.num_colors + i, endpoint=True)
        for thd, widget in zip(thds[i:self.num_colors][::-1], self.values_widgets):
            # Remove bounds, set values, update bounds
            widget.start = None
            widget.end = None
            widget.value = thd
            self._update_bounds()

    @param.depends('mode', watch=True)
    def _mode_updated(self):
        if self.mode == 'Discrete':
            self.box_insert_after('num_colors', 'otsu_thd')
            #self.otsu_thd.constant = False
        elif self.mode == 'Continuous':
            self.box_pop('otsu_thd')
        elif self.mode == 'ColorMap':
            self.num_colors = 2
            #todo adjust add/ remove color widgets methods
        self.param.trigger('num_colors')

    @param.depends('values', 'colors', 'target', 'quantity', watch=True)
    def _get_colors(self):
        # todo or?
        if np.all(self.values == 0):
            return
        elif np.any(np.diff(self.values) > 0):  # Skip applying colors when not strictly monotonic descending
            return
        elif not self.target:
            return

        y_vals = self.parent.sources[self.target][self.quantity]  # full array including nan entries

        if self.num_colors == 1:
            colors = np.full(len(y_vals), fill_value=self.colors[0], dtype='U7')
            colors[np.isnan(y_vals)] = np.nan
        elif self.mode == 'Discrete':
            full_thds = [-np.inf] + self.values[::-1] + [np.inf]
            colors = np.full(len(y_vals), fill_value=np.nan, dtype='U7')
            for lower, upper, color in zip(full_thds[:-1], full_thds[1:], self.colors[::-1]):
                b = (y_vals > lower) & (y_vals <= upper)
                colors[b] = color
        elif self.mode == 'Continuous':
            func = np.log if self.log_space else lambda x: x
            vals_space = (func(self.values))  # values in log space depending on setting
            norm = plt.Normalize(vals_space[-1], vals_space[0], clip=True)
            nodes = norm(vals_space[::-1])
            cmap = mpl.colors.LinearSegmentedColormap.from_list("custom_cmap", list(zip(nodes, self.colors[::-1])))

            try:
                colors_rgba = cmap(norm(func(y_vals)))
                colors = np.array([rgb_to_hex(int(r*255), int(g*255), int(b*255)) for r, g, b, a in colors_rgba])
                colors[np.isnan(y_vals)] = np.nan
            except ValueError as err:
                self.parent.logger.warning(err)
                return

        self.parent.sources[self.target].source.data['color'] = colors  # this triggers an update of the graph

    @param.depends('num_colors', watch=True)
    def _update_num_colors(self):
        while len(self.colors_widgets) != self.num_colors:
            if len(self.colors_widgets) > self.num_colors:
                self._remove_color()
            elif len(self.colors_widgets) < self.num_colors:
                self._add_color()
        self.param.trigger('colors')

    @param.depends('num_colors', watch=True)
    def _update_num_values(self):
        diff = 1 if self.mode == 'Discrete' else 0
        while len(self.values_widgets) != self.num_colors - diff:
            if len(self.values_widgets) > self.num_colors - diff:
                self._remove_value()
            elif len(self.values_widgets) < self.num_colors - diff:
                self._add_value()

        self._update_bounds()
        self.param.trigger('values')

    def _add_value(self):
        try:
            first_value = self.values_widgets[-1].value
        except IndexError:
            first_value = 0

        default = float(first_value - 1)
        self.values.append(default)

        name = 'Threshold {}'.format(len(self.values_widgets) + 1)
        widget = pn.widgets.FloatInput(name=name, value=default)
        self.values_widgets.append(widget)
        i = len(self.values_widgets) + self.box_index('show_thds')
        self._box.insert(i, widget)
        widget.param.watch(self._value_event, ['value'])

    def _remove_value(self):
        widget = self.values_widgets.pop(-1)
        self.box_pop(widget)
        self.values.pop()

        [widget.param.unwatch(watcher) for watcher in widget.param._watchers]
        del widget

    def _add_color(self):
        try:
            default = DEFAULT_CLASS_COLORS[len(self.colors_widgets)]
        except IndexError:
            default = "#"+''.join(np.random.choice(list('0123456789abcdef'), 6))

        self.colors.append(default)
        widget = pn.widgets.ColorPicker(value=default)
        self.colors_widgets.append(widget)
        i = len(self.values_widgets) + len(self.colors_widgets) + self.box_index('show_thds')
        self._box.insert(i, widget)
        widget.param.watch(self._color_event, ['value'])

    def _remove_color(self):
        widget = self.colors_widgets.pop(-1)
        self.colors.pop()
        self.box_pop(widget)
        [widget.param.unwatch(watcher) for watcher in widget.param._watchers]
        del widget

    def _color_event(self, *events):
        for event in events:
            idx = list(self.colors_widgets).index(event.obj)
            self.colors[idx] = event.new


        #todo callback?
        self.param.trigger('colors')

    def _value_event(self, *events):
        """triggers when a single value gets changed"""
        for event in events:
            idx = list(self.values_widgets).index(event.obj)
            self.values[idx] = event.new

        self._update_bounds()
        self.param.trigger('values')

    def _update_bounds(self):
        for i, widget in enumerate(self.values_widgets):
            if i > 0:
                prev_value = float(self.values_widgets[i - 1].value)
                widget.end = np.nextafter(prev_value, prev_value - 1)
            else:
                widget.end = None

            if i < len(self.values_widgets) - 1:
                next_value = float(self.values_widgets[i + 1].value)
                widget.start = np.nextafter(next_value, next_value + 1)
            else:
                widget.start = None


class FileExportControl(ControlPanel):
    # todo check if docstring is true
    """
    This controller allows users to export and download datasets.

    All datasets can be exported as .txt tables.
    'Mappable' datasets (with r_number column) can be exported as .pml pymol script, which colors protein structures
    based on their 'color' column.

    """

    header = "File Export"
    target = param.Selector(label='Target dataset', doc='Name of the dataset to export')
    #todo add color param an dlink with protein viewer color

    def __init__(self, parent, **param):
        self.export_linear_download = pn.widgets.FileDownload(filename='<no data>', callback=self.linear_export_callback)
        self.pml_script_download = pn.widgets.FileDownload(filename='<no data>', callback=self.pml_export_callback)
        super(FileExportControl, self).__init__(parent, **param)

        self.parent.param.watch(self._sources_updated, ['sources'])
        try:  # todo write function that does the try/excepting + warnings (or also mixins?)
            self.parent.param.watch(self._series_updated, ['series'])
        except ValueError:
            pass

    def make_list(self):
        self._widget_dict.update(export_linear_download=self.export_linear_download, pml_script_download=self.pml_script_download)
        return super(FileExportControl, self).make_list()

    def _sources_updated(self, *events):
        objects = list(self.parent.sources.keys())
        self.param['target'].objects = objects

        if not self.target and objects:
            self.target = objects[0]

    def _series_updated(self, *events):
        self.c_term = int(self.parent.series.cov.protein.c_term)

    def _make_pml(self, target):
        try:
            #todo add no coverage field and link to the other no coverage field
            no_coverage = self.parent.control_panels['ProteinViewControl'].no_coverage
        except KeyError:
            no_coverage = '#8c8c8c'
            self.parent.logger.warning('No coverage color found, using default grey')

        try:
            script = colors_to_pymol(self.export_dict['r_number'], self.export_dict['color'],
                                     c_term=self.parent.series.c_term, no_coverage=no_coverage)
            return script
        except KeyError:
            return None

    @property
    def export_dict(self):
        return self.export_data_source.source.data

    @property
    def export_data_source(self):
        return self.parent.sources[self.target]

    @pn.depends('target', watch=True)
    def _update_filename(self):
        #todo subclass and split
        self.export_linear_download.filename = self.parent.series.state + '_' + self.target + '_linear.txt'
        if 'mapping' in self.export_data_source.tags:
            self.pml_script_download.filename = self.parent.series.state + '_' + self.target + '_pymol.pml'
            # self.pml_script_download.disabled = False
        else:
            self.pml_script_download.filename = 'Not Available'
            # self.pml_script_download.disabled = True # Enable/disable currently bugged:

    @pn.depends('target')
    def pml_export_callback(self):
        if self.target:
            io = StringIO()
            io.write('# ' + VERSION_STRING + ' \n')
            script = self._make_pml(self.target)
            try:
                io.write(script)
                io.seek(0)
                return io
            except TypeError:
                return None
        else:
            return None

    @pn.depends('target')  # param.depends?
    def linear_export_callback(self):
        io = StringIO()
        io.write('# ' + VERSION_STRING + ' \n')

        if self.target:
            export_dict = {k: np.array(v) for k, v in self.parent.sources[self.target].source.data.items()}  #todo generalize export
            dtype = [(name, arr.dtype) for name, arr in export_dict.items()]
            export_data = np.empty_like(next(iter(export_dict.values())), dtype=dtype)
            for name, arr in export_dict.items():
                export_data[name] = arr

            fmt, header = fmt_export(export_data)
            np.savetxt(io, export_data, fmt=fmt, header=header)

            io.seek(0)
            return io
        else:
            return None


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



class ProteinViewControl(ControlPanel):
    """
    This controller allows users control the Protein view figure.
    Structures can be specified either by RCSB ID or uploading a .pdb file.

    Colors are assigned according to 'color' column of the selected dataset.
    """

    header = 'Protein Viewer'
    accepted_tags = ['mapping']

    target_dataset = param.Selector(doc='Name of the dataset to apply coloring from')
    input_option = param.Selector(default='Upload File', objects=['Upload File', 'RCSB PDB'],
                                  doc='Choose wheter to upload .pdb file or directly download from RCSB PDB.')
    rcsb_id = param.String(doc='RCSB PDB identifier of protein entry to download and visualize.')
    #load_structure = param.Action(lambda self: self._load_structure())
    no_coverage = param.Color(default='#8c8c8c', doc='Color to use for regions of no coverage.')
    representation = param.Selector(default='cartoon',
                                    objects=['backbone', 'ball+stick', 'cartoon', 'hyperball', 'licorice',
                                             'ribbon', 'rope', 'spacefill', 'surface'],
                                    doc='Representation to use to render the protein.')
    spin = param.Boolean(default=False, doc='Rotate the protein around an axis.')

    def __init__(self, parent, **params):
        self.file_input = pn.widgets.FileInput(accept='.pdb')
        super(ProteinViewControl, self).__init__(parent, **params)

        self.parent.param.watch(self._parent_sources_updated, ['sources'])
        self.input_option = 'RCSB PDB'

    def make_list(self):
        lst = super().make_list()
        lst.pop(2)  # Remove RCSB ID input field?
        lst.insert(2, self.file_input)  # add File input widget
        return lst

    def _parent_sources_updated(self, *events):
        #todo  this line repeats, put in base class
        data_sources = [k for k, src in self.parent.sources.items() if src.resolve_tags(self.accepted_tags)]
        self.param['target_dataset'].objects = data_sources
        if not self.target_dataset and data_sources:
            self.target_dataset = data_sources[0]

    @param.depends('input_option', watch=True)
    def _update_input_option(self):
        if self.input_option == 'Upload File':
            self.box_pop('rcsb_id')
            self.box_insert_after('input_option', self.file_input)
        elif self.input_option == 'RCSB PDB':
            self.box_pop(self.file_input)
            self.box_insert_after('input_option', 'rcsb_id')

        elif self.input_option == 'RCSB PDB':
            self.ngl_html.rcsb_id = self.rcsb_id


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

    def __init__(self, parent, **params):
        super(DeveloperControl, self).__init__(parent, **params)

    def _action_test_logging(self):
        self.parent.logger.debug('TEST DEBUG MESSAGE')
        for i in range(20):
            self.parent.logger.info('dit is een test123')

    def _action_print(self):
        print(self.parent.doc)

    def _action_break(self):
        main_ctrl = self.parent
        control_panels = main_ctrl.control_panels
        figure_panels = main_ctrl.figure_panels
        sources = main_ctrl.sources

        print('Time for a break')

    def _action_test(self):
        from pathlib import Path
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