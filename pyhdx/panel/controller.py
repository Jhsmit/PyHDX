
from .log import setup_custom_logger
from .base import ControlPanel
from .panels import FileInputPanel, CoveragePanel, RateConstantPanel

from .fig_panels import CoverageFigure, RateFigure
from pyhdx.pyhdx import PeptideCSVFile, KineticsSeries
from pyhdx.fitting import KineticsFitting
from pyhdx.fileIO import read_dynamx
from pyhdx.support import get_constant_blocks, get_reduced_blocks

logger = setup_custom_logger('root')
logger.debug('main message')


import param
import panel as pn
from jinja2 import Environment, FileSystemLoader
import holoviews as hv
import os
import numpy as np
from numpy.lib.recfunctions import stack_arrays
from io import StringIO

import matplotlib
matplotlib.use('agg') # for panel mpl support

from bokeh.util.serialization import make_globally_unique_id


pth = os.path.dirname(__file__)

env = Environment(loader=FileSystemLoader(pth))


class Controller(param.Parameterized):
    """
    controller for main panels layout
    and has panels for each tabin the main layout

    """

    data = param.Array()  # might not be needed, in favour of peptides
    rates = param.Array(doc='Output rates data')
    peptides = param.ClassSelector(PeptideCSVFile)  #class with all peptides to be considered
    series = param.ClassSelector(KineticsSeries)
    fitting = param.ClassSelector(KineticsFitting)

    def __init__(self, template, panels, **params):
        super(Controller, self).__init__(**params)
        template = env.get_template('template.html')

        tmpl = pn.Template(template=template)
     #   tmpl.nb_template.globals['get_id'] = make_globally_unique_id


        # Controllers
        self.fileinput = FileInputControl(self)
        self.coverage = CoverageControl(self)#CoveragePanel(self)
        self.fit_control = FittingControl(self)
        self.rate_panel = RateConstantPanel(self)


        #Figures
        self.coverage_figure = CoverageFigure(self, [self.coverage])  #parent, [controllers]
        self.rate_figure = RateFigure(self, [self.fit_control]) # parent, [controllers]

        # tmpl = pn.Template(template)
        tmpl.add_panel('input', self.fileinput.panel)
        tmpl.add_panel('coverage', self.coverage.panel)
        tmpl.add_panel('fitting', self.fit_control.panel)
        tmpl.add_panel('scene3d', self.coverage_figure.panel)
        tmpl.add_panel('slice_j', self.rate_figure.panel)
        tmpl.add_panel('slice_k', hv.Curve([1, 2, 3]))


        #tmpl.add_panel('B', hv.Curve([1, 2, 3]))

        self.template = tmpl
      #  self.panels = [panel(self) for panel in panels]

    @param.depends('series', watch=True)
    def _series_changed(self):
        # This is triggered if the fileinput child panel yields a new KineticSeries
        print('series changed')

        self.fitting = KineticsFitting(self.series)
        #todo add errors here
        rate_fields = ['fit1', 'fit1_r1', 'fit1_r2', 'fit2', 'fit2_r1', 'fit2_r2']
        rates = np.zeros(self.series.cov.prot_len,
                              dtype=[('r_number', int)] + [(name, float ) for name in rate_fields])
        rates['r_number'] = self.series.cov.r_number

        self.rates = rates  # this assignement triggers downstream watchers

    @property
    def servable(self):
        return self.template.servable

    @param.depends('data')
    def _test(self):
        print("hoi, data changed")


class FileInputControl(ControlPanel):
    add_button = param.Action(lambda self: self._action_add(), doc='Add File', label='Add File')
    clear_button = param.Action(lambda self: self._action_clear(), doc='Clear files', label='Clear Files')
    drop_first = param.Integer(1, bounds=(0, None))
    load_button = param.Action(lambda self: self._action_load(), doc='Load Files', label='Load Files')

    control_state = param.Selector(doc='State for the control condition', label='Control State')
    control_exposure = param.Selector(doc='Exposure for control condition', label='Control exposure')

    exp_state = param.Selector(doc='State for selected experiment', label='Experiment State')
    exp_times = param.ListSelector(default=[], objects=[''])

    parse_button = param.Action(lambda self: self._action_parse(), doc='Parse', label='Parse')

    def __init__(self, parent, **params):
        super(FileInputControl, self).__init__(parent, **params)
        self.file_selectors_column = pn.Column(*[pn.widgets.FileInput(accept='.csv')])

    def _action_add(self):
        print('action_add')
        widget = pn.widgets.FileInput(accept='.csv')
        self.file_selectors_column.append(widget)

    def _action_clear(self):
        print('action clear')
        self.file_selectors_column.clear()
        self._action_add()

    def _action_load(self):
        print('action load')
        data_list = []
        for file_selector in self.file_selectors_column:
            if file_selector.value is not None:
                s_io = StringIO(file_selector.value.decode('UTF-8'))
                data = read_dynamx(s_io)
                data_list.append(data)

        combined = stack_arrays(data_list, asrecarray=True, usemask=False, autoconvert=True)

        self.parent.data = combined
        self.parent.peptides = PeptideCSVFile(self.parent.data, drop_first=self.drop_first)

        states = list(np.unique(self.parent.peptides.data['state']))
        self.param['control_state'].objects = states
        self.control_state = states[0]

    def _action_parse(self):
        print('parse action')
        self.parent.peptides.set_control((self.control_state, self.control_exposure), remove_nan=True)
        data_states = self.parent.peptides.data[self.parent.peptides.data['state'] == self.exp_state]
        data = data_states[np.isin(data_states['exposure'], self.exp_times)]
        series = KineticsSeries(data)
        series.make_uniform()  #TODO add gui control for this

        self.parent.series = series

    @param.depends('control_state', watch=True)
    def _update_control_exposure(self):
        b = self.parent.peptides.data['state'] == self.control_state
        data = self.parent.peptides.data[b]
        exposures = list(np.unique(data['exposure']))
        self.param['control_exposure'].objects = exposures
        self.control_exposure = exposures[0]

    @param.depends('control_state', 'control_exposure', watch=True)
    def _update_experiment(self):
        # r = str(np.random.rand())
        # self.param['exp_state'].objects = [r]
        # self.exp_state = r

        pm_dict = self.parent.peptides.return_by_name(self.control_state, self.control_exposure)
        states = list(np.unique([v.state for v in pm_dict.values()]))
        self.param['exp_state'].objects = states
        self.exp_state = states[0]

    @param.depends('exp_state', watch=True)
    def _update_experiment_exposure(self):
        b = self.parent.data['state'] == self.exp_state
        exposures = list(np.unique(self.parent.data['exposure'][b]))
        exposures.sort()
        self.param['exp_times'].objects = exposures
        self.exp_times = exposures

    @property
    def panel(self):
        params = pn.panel(self.param)
        col = pn.Column(*[self.file_selectors_column, *params[1:]])

        return pn.WidgetBox(col, pn.layout.VSpacer(), css_classes=['widget-box', 'custom-wbox'], sizing_mode='stretch_height')


class CoverageControl(ControlPanel):
    wrap = param.Integer(25, bounds=(0, None), doc='Number of peptides vertically before moving to the next row') # todo auto?
    aa_per_subplot = param.Integer(100, label='Amino acids per subplot')
    update = param.Action(lambda self: self.param.trigger('update'), label='Update')  #Triggers redraw of the figure (might not need this in favour of directly triggering from wrap and aa per subplot
    labels = param.Boolean(False, label='Labels')
    index = param.Integer(0, bounds=(0, None), doc='Current index of coverage plot in time')

    def __init__(self, parent, **params):
        super(CoverageControl, self).__init__(parent, **params)
        self.exposure_str = pn.widgets.StaticText(name='Exposure', value='0') # todo update to some param?
        self.parent.param.watch(self._update_series, ['series'])

    @property
    def peptide_measurement(self):
        return self.parent.series[self.index]

    def _update_series(self, event):
        print('coverage new series update index bounds')
        self.param['index'].bounds = (0, len(event.new) - 1)
        self.exposure_str.value = str(self.peptide_measurement.exposure)

    @param.depends('index', watch=True)
    def _update_index(self):
        self.exposure_str.value = str(self.peptide_measurement.exposure)

    @property
    def panel(self):
        col = pn.Column(self.param, self.exposure_str)
        #p = pn.Param(self.param, widgets={'file': pn.widgets.FileInput}) for exposure
        return pn.WidgetBox(col, pn.layout.VSpacer(), css_classes=['widget-box', 'custom-wbox'],
                            sizing_mode='stretch_height')


class FittingControl(ControlPanel):
    chisq_thd = param.Number(20, doc='Threshold for chi2 to switch to Differential evolution')
    r_max = param.Number(27, doc='Ceil value for rates')  # Update this value

    do_fit1 = param.Action(lambda self: self._action_fit1())
    block_mode = param.ObjectSelector(default='reduced', objects=['reduced', 'original', 'constant'])

    #todo generate from func signature?
    #block mode reduced params
    max_combine = param.Integer(2, doc='Neighbouring blocks up to and including this size are merged together')
    max_join = param.Integer(5, doc='Blocks up to and including this size are joined with their smallest neighbour')

    #constant block params
    block_size = param.Integer(10, doc='Size of the blocks in constant blocks mode')
    initial_block = param.Integer(5, doc='Size of the initial block in constant block mode')
    show_blocks = param.Boolean(False, doc='Show bounds of blocks in graph')

    do_fit2 = param.Action(lambda self: self._action_fit2(), constant=True)

    def __init__(self, parent, **params):
        super(FittingControl, self).__init__(parent, **params)

        self.block_column = pn.Column(*[self.param[key] for key in ['max_combine', 'max_join']])

        self.parent.param.watch(self._update_series, ['series'])

    def _update_series(self, *events):
        self.r_max = np.log(1 - 0.98) / -self.parent.series.times[1]  # todo user input 0.98

    def _action_fit1(self):
        print('fitting 1')
        #todo context manager?
        self.param['do_fit1'].constant = True
        self.param['do_fit2'].constant = True

        result = self.parent.fitting.weighted_avg_fit(chisq_thd=self.chisq_thd)
        rates = np.vstack([1 / result.get_param('tau1'), 1 / result.get_param('tau2')])

        names = ['fit1', 'fit1_r1', 'fit1_r2']
        rates_out = np.empty_like(self.parent.rates,
                                  dtype=[(name, float) for name in names])
        rates_out['fit1'] = result.rate
        rates_out['fit1_r1'] = rates.min(axis=0)
        rates_out['fit1_r2'] = rates.max(axis=0)

        self.parent.rates[names] = rates_out[names]  # assigning fields doesnt seem to trigger replot
        self.parent.param.trigger('rates')
        #self._renew(None)  #manual trigger

        self.param['do_fit1'].constant = False
        self.param['do_fit2'].constant = False

    def _action_fit2(self):
        print('fitting 2')
        #todo context manager
        self.param['do_fit1'].constant = True
        self.param['do_fit2'].constant = True

        r_number, fit2_rate = self.parent.fitting.lsq_fit_blocks(self.parent.rates, **self.fit_kwargs)
        self.parent.rates['fit2'] = fit2_rate
        self.parent.param.trigger('rates')

  #      self._renew(None)  # manual trigger

        self.param['do_fit1'].constant = False
        self.param['do_fit2'].constant = False

    @property
    def fit_kwargs(self):
        if self.block_mode == 'reduced':
            fit_kwargs = {'block_func': get_reduced_blocks, 'max_combine': self.max_combine, 'max_join': self.max_join}
        elif self.block_mode == 'original':
            fit_kwargs = {'block_func': lambda series, **kwargs: series.cov.block_length}
        elif self.block_mode == 'constant':
            fit_kwargs = {'block_func': get_constant_blocks, 'block_size': self.block_size, 'initial_block': self.initial_block}
        return fit_kwargs

    @param.depends('block_mode', watch=True)
    def _update_block_mode(self):
        print('block mode updated')
        if self.block_mode == 'reduced':
            self.block_column.clear()
            [self.block_column.append(self.param[key]) for key in ['max_combine', 'max_join']]
        elif self.block_mode == 'original':
            self.block_column.clear()
        elif self.block_mode == 'constant':
            self.block_column.clear()
            [self.block_column.append(self.param[key]) for key in ['block_size', 'initial_block']]

    @property
    def panel(self):
        par1 = ['chisq_thd', 'r_max', 'do_fit1', 'block_mode']
        par2 = ['do_fit2']

        pars = [self.param[key] for key in par1] + [self.block_column] + [self.param[key] for key in par2]
        return pn.WidgetBox(*pars)
