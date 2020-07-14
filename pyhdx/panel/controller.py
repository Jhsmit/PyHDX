from .log import setup_custom_logger
from .base import ControlPanel, DEFAULT_COLORS, DEFAULT_CLASS_COLORS
from .fig_panels import CoverageFigure, RateFigure, ProteinFigure, FitResultFigure, PFactFigure
from pyhdx.pyhdx import PeptideMasterTable, KineticsSeries
from pyhdx.fitting import KineticsFitting
from pyhdx.fileIO import read_dynamx
from pyhdx.support import get_constant_blocks, get_reduced_blocks, get_original_blocks, fmt_export, np_from_txt, \
    autowrap, colors_to_pymol

from pyhdx import VERSION_STRING

logger = setup_custom_logger('root')
logger.debug('main message')

from scipy import constants
import param
import panel as pn
from jinja2 import Environment, FileSystemLoader
#import holoviews as hv  #todo remove dependency
import os
import numpy as np
from skimage.filters import threshold_multiotsu
from numpy.lib.recfunctions import stack_arrays, append_fields
from .components import ASyncProgressBar
from io import StringIO, BytesIO
from tornado.ioloop import IOLoop
from bokeh.plotting import curdoc
from bokeh.document import without_document_lock
import matplotlib
matplotlib.use('agg') # for panel mpl support
from functools import partial
#from .widgets import NumericInput
from bokeh.models import ColumnDataSource

#dev only
import pickle


from bokeh.util.serialization import make_globally_unique_id
pth = os.path.dirname(__file__)

env = Environment(loader=FileSystemLoader(pth))

# todo dict comprehension

#refactor rates to columndatasource?
dic = {'rates': np.zeros(0, dtype=[('r_number', int), ('rate', float)]),
       'fitresult': None}

empty_results = {
    'fit1': dic.copy(),
    'fit2': dic.copy()
}


class Controller(param.Parameterized):
    """
    controller for main panels layout
    and has panels for each tabin the main layout

    """

    data = param.Array()  # might not be needed, in favour of peptides
    #rates = param.Array(doc='Output rates data')
    fit_results = param.Dict(empty_results)
    sources = param.Dict({}, doc='Dictionary of ColumnDataSources available for plotting')
    rate_colors = param.Dict({})  #probably not used
    peptides = param.ClassSelector(PeptideMasterTable, doc='Master list of all peptides')
    series = param.ClassSelector(KineticsSeries, doc='Currently selected kinetic series of peptides')
    fitting = param.ClassSelector(KineticsFitting)

    def __init__(self, template, panels, cluster=None, **params):
        super(Controller, self).__init__(**params)
        template = env.get_template('template.html')
        self.cluster = cluster
        self.doc = pn.state.curdoc
        tmpl = pn.Template(template=template)
     #   tmpl.nb_template.globals['get_id'] = make_globally_unique_id

        # Controllers
        self.fileinput = FileInputControl(self)
        self.coverage = CoverageControl(self)#CoveragePanel(self)
        self.fit_control = FittingControl(self)
        self.tf_fit_control = TFFitControl(self)
        self.fit_quality = FittingQuality(self)
        #self.rate_panel = RateConstantPanel(self)
        self.classification_panel = ClassificationControl(self)
        self.file_export = FileExportPanel(self)
        self.options = OptionsPanel(self)
        self.dev = DeveloperPanel(self)

        #Figures
        self.coverage_figure = CoverageFigure(self, [self.coverage, self.fit_control])  #parent, [controllers]
        self.rate_figure = RateFigure(self, [self.fit_control, self.classification_panel]) # parent, [controllers]  #todo parse as kwargs
        self.pfact_figure = PFactFigure(self, [self.fit_control, self.classification_panel])

        self.fit_result_figure = FitResultFigure(self, [self.fit_quality])
        self.protein_figure = ProteinFigure(self, [])


        #setup options  #todo automate figure out cross dependencies
        self.options.cov_fig_panel = self.coverage_figure
        self.options.rate_fig_panel = self.rate_figure
        self.options.coverage_ctrl = self.coverage

        # tmpl = pn.Template(template)
        tmpl.add_panel('input', self.fileinput.panel)
        tmpl.add_panel('coverage', self.coverage.panel)
        tmpl.add_panel('fitting', self.fit_control.panel)
        tmpl.add_panel('tf_fit', self.tf_fit_control.panel)
        tmpl.add_panel('fit_quality', self.fit_quality.panel)
        tmpl.add_panel('classification', self.classification_panel.panel)
        tmpl.add_panel('file_export', self.file_export.panel)
        tmpl.add_panel('options', self.options.panel)
        tmpl.add_panel('dev', self.dev.panel)

        tmpl.add_panel('coverage_fig', self.coverage_figure.panel)

        tmpl.add_panel('rate_fig', self.rate_figure.panel)
        tmpl.add_panel('pfact_fig', self.pfact_figure.panel)
        tmpl.add_panel('fitres_fig', self.fit_result_figure.panel)
        tmpl.add_panel('slice_k', self.protein_figure.panel)
        #tmpl.add_panel('B', hv.Curve([1, 2, 3]))

        self.app = tmpl
      #  self.panels = [panel(self) for panel in panels]

    @param.depends('series', watch=True)
    def _series_changed(self):
        # This is triggered if the fileinput child panel yields a new KineticSeries
        print('series changed')

        self.fitting = KineticsFitting(self.series, cluster=self.cluster)
        for key in ['fit1', 'fit2']:    # todo this list of names somewhere?
            self.rate_colors[key] = [DEFAULT_COLORS[key]]*len(self.series.cov.r_number)
        self.param.trigger('rate_colors')

        # #todo add errors here
        # rate_fields = ['fit1', 'fit1_r1', 'fit1_r2', 'fit2', 'fit2_r1', 'fit2_r2']
        # color_fields = ['fit1_color', 'fit2_color']
        # dtype = [('r_number', int)] + [(name, float) for name in rate_fields] + [(name, 'U7') for name in color_fields]
        # rates = np.zeros(self.series.cov.prot_len, dtype=dtype)
        # rates['r_number'] = self.series.cov.r_number
        # rates['fit1_color'][:] = 'blue'
        # rates['fit2_color'][:] = 'red'
        #
        # self.rates = rates  # this assignement triggers downstream watchers? manual trigger?

    def servable(self):

        js_files = {'jquery': 'https://code.jquery.com/jquery-1.11.1.min.js',
                    'goldenlayout': 'https://golden-layout.com/files/latest/js/goldenlayout.min.js',
                    'ngl': 'https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.33/dist/ngl.js'}
        css_files = ['https://golden-layout.com/files/latest/css/goldenlayout-base.css',
                     'https://golden-layout.com/files/latest/css/goldenlayout-dark-theme.css']

        css = '''
        .custom-wbox > div.bk {
            padding-right: 10px;
        }
        .scrollable {
            overflow: auto !important;
        }
        '''
        pn.extension(js_files=js_files, raw_css=[css], css_files=css_files)

        return self.app.servable()

    def serve(self, **kwargs):
        js_files = {'jquery': 'https://code.jquery.com/jquery-1.11.1.min.js',
                    'goldenlayout': 'https://golden-layout.com/files/latest/js/goldenlayout.min.js',
                    'ngl': 'https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.33/dist/ngl.js'}
        css_files = ['https://golden-layout.com/files/latest/css/goldenlayout-base.css',
                     'https://golden-layout.com/files/latest/css/goldenlayout-dark-theme.css']

        css = '''
        .custom-wbox > div.bk {
            padding-right: 10px;
        }
        .scrollable {
            overflow: auto !important;
        }
        '''
        pn.extension(js_files=js_files, raw_css=[css], css_files=css_files)
        pn.serve(self.app, **kwargs)

    def get_rate_file_export(self):
        fmt, header = fmt_export(self.rates)
        s = StringIO()
        np.savetxt(s, fmt=fmt, header=header)

    @param.depends('data')
    def _test(self):
        print("hoi, data changed")


class FileInputControl(ControlPanel):
    header = 'Input'

    add_button = param.Action(lambda self: self._action_add(), doc='Add File', label='Add File')
    clear_button = param.Action(lambda self: self._action_clear(), doc='Clear files', label='Clear Files')
    drop_first = param.Integer(1, bounds=(0, None))
    ignore_prolines = param.Boolean(True, doc='Set to True to ignore prolines in the sequence')
    load_button = param.Action(lambda self: self._action_load(), doc='Load Files', label='Load Files')

    norm_mode = param.Selector(doc='Select method of normalization', label='Norm mode', objects=['Exp', 'Theory'])

    norm_state = param.Selector(doc='State used to normalize uptake', label='Norm State')
    norm_exposure = param.Selector(doc='Exposure used to normalize uptake', label='Norm exposure')
    be_percent = param.Number(28., bounds=(0, 100), doc='Percentage of exchangeable deuteriums which backexchange',
                              label='Back exchange percentage')

    zero_state = param.Selector(doc='State used to zero uptake', label='Zero state')
    zero_exposure = param.Selector(doc='Exposure used to zero uptake', label='Zero exposure')

    exp_state = param.Selector(doc='State for selected experiment', label='Experiment State')
    exp_exposures = param.ListSelector(default=[], objects=[''], label='Experiment Exposures')

    parse_button = param.Action(lambda self: self._action_parse(), doc='Parse', label='Parse')

    def __init__(self, parent, **params):
        self.file_selectors = [pn.widgets.FileInput(accept='.csv')]
        super(FileInputControl, self).__init__(parent, **params)

    def make_dict(self):
        return self.generate_widgets(norm_mode=pn.widgets.RadioButtonGroup, be_percent=pn.widgets.LiteralInput)

    def make_list(self):
        parameters = ['add_button', 'clear_button', 'drop_first', 'ignore_prolines', 'load_button',
                      'norm_mode', 'norm_state', 'norm_exposure',  'zero_state', 'zero_exposure', 'exp_state',
                      'exp_exposures', 'parse_button']
        first_widgets = list([self._widget_dict[par] for par in parameters])
        return self.file_selectors + first_widgets

    def _action_add(self):
        print('action_add')
        widget = pn.widgets.FileInput(accept='.csv')
        i = len(self.file_selectors) + 1 # position to insert the new file selector into the widget box
        self.file_selectors.append(widget)
        self._box.insert(i, widget)

    def _action_clear(self):
        print('action clear')

        while self.file_selectors:
            fs = self.file_selectors.pop()
            #todo allow popping/locking with both widgets and parameter names?
            idx = list(self._box).index(fs)
            self._box.pop(idx)
        self._action_add()

    def _action_load(self):
        print('action load')
        data_list = []
        for file_selector in self.file_selectors:
            if file_selector.value is not None:
                s_io = StringIO(file_selector.value.decode('UTF-8'))
                data = read_dynamx(s_io)
                data_list.append(data)

        combined = stack_arrays(data_list, asrecarray=True, usemask=False, autoconvert=True)

        self.parent.data = combined
        self.parent.peptides = PeptideMasterTable(self.parent.data,
                                                  drop_first=self.drop_first, ignore_prolines=self.ignore_prolines)

        states = list(np.unique(self.parent.peptides.data['state']))
        self.param['norm_state'].objects = states
        self.norm_state = states[0]
        self.param['zero_state'].objects = ['None'] + states
        self.zero_state = 'None'

    def _action_parse(self):
        print('parse action')
        if self.norm_mode == 'Exp':
            control_0 = (self.zero_state, self.zero_exposure) if self.zero_state != 'None' else None
            self.parent.peptides.set_control((self.norm_state, self.norm_exposure), control_0=control_0)
        elif self.norm_mode == 'Theory':
            self.parent.peptides.set_backexchange(self.be_percent)

        data_states = self.parent.peptides.data[self.parent.peptides.data['state'] == self.exp_state]
        data = data_states[np.isin(data_states['exposure'], self.exp_exposures)]

        # states = self.parent.peptides.groupby_state()
        # series = states[self.exp_state]
        # series.make_uniform()

        # b = np.isin(series.full_data['exposure'], self.exp_exposures)
        # data = series.full_data[b].copy()

        #series = KineticsSeries(data)
        #series.make_uniform()  #TODO add gui control for this

        series = KineticsSeries(data)
        series.make_uniform()
        self.parent.series = series

    @param.depends('norm_mode', watch=True)
    def _update_norm_mode(self):

        if self.norm_mode == 'Exp':
            self.box_pop('be_percent')
            self.box_insert_after('norm_mode', 'norm_state')
            self.box_insert_after('norm_state', 'norm_exposure')
            self.box_insert_after('norm_exposure', 'zero_state')
            self.box_insert_after('zero_state', 'zero_exposure')

            #self._update_experiment()  dont think this is needed
        elif self.norm_mode == 'Theory':
            self.box_pop('norm_state')
            self.box_pop('norm_exposure')
            self.box_pop('zero_state')
            self.box_pop('zero_exposure')
            self.box_insert_after('norm_mode', 'be_percent')

            try:
                states = np.unique(self.parent.data['state'])
                self.param['exp_state'].objects = states
                self.exp_state = states[0] if not self.exp_state else self.exp_state
            except TypeError:
                pass

    @param.depends('norm_state', watch=True)
    def _update_norm_exposure(self):
        b = self.parent.peptides.data['state'] == self.norm_state
        data = self.parent.peptides.data[b]
        exposures = list(np.unique(data['exposure']))
        self.param['norm_exposure'].objects = exposures
        if exposures:
            self.norm_exposure = exposures[0]

    @param.depends('zero_state', watch=True)
    def _update_zero_exposure(self):
        b = self.parent.peptides.data['state'] == self.zero_state
        data = self.parent.peptides.data[b]
        exposures = list(np.unique(data['exposure']))
        self.param['zero_exposure'].objects = exposures
        if exposures:
            self.control_exposure = exposures[0]

    @param.depends('norm_state', 'norm_exposure', watch=True)
    def _update_experiment(self):
        # r = str(np.random.rand())
        # self.param['exp_state'].objects = [r]
        # self.exp_state = r
        #TODO THIS needs to be updated to also incorporate the zero
        print(self.norm_state, self.norm_exposure)
        pm_dict = self.parent.peptides.return_by_name(self.norm_state, self.norm_exposure)
        states = list(np.unique([v.state for v in pm_dict.values()]))
        self.param['exp_state'].objects = states
        self.exp_state = states[0] if not self.exp_state else self.exp_state

    @param.depends('exp_state', watch=True)
    def _update_experiment_exposure(self):
        b = self.parent.data['state'] == self.exp_state
        exposures = list(np.unique(self.parent.data['exposure'][b]))
        exposures.sort()
        self.param['exp_exposures'].objects = exposures  #todo refactor exposures
        self.exp_exposures = exposures


class CoverageControl(ControlPanel):
    header = 'Coverage'

    wrap = param.Integer(25, bounds=(0, None), doc='Number of peptides vertically before moving to the next row') # todo auto?
    aa_per_subplot = param.Integer(100, label='Amino acids per subplot')
    labels = param.Boolean(False, label='Labels')
    index = param.Integer(0, bounds=(0, 10), doc='Current index of coverage plot in time')

    def __init__(self, parent, **params):
        self.exposure_str = pn.widgets.StaticText(name='Exposure', value='0') # todo update to some param?
        super(CoverageControl, self).__init__(parent, **params)
        self.parent.param.watch(self._update_series, ['series'])

    def make_list(self):
        lst = super(CoverageControl, self).make_list()
        return lst + [self.exposure_str]

    def make_dict(self):
        return self.generate_widgets(index=pn.widgets.IntSlider)

    @property
    def peptide_measurement(self):
        if self.parent.series is not None:
            return self.parent.series[self.index]
        else:
            return None

    def _update_series(self, event):
        print('coverage new series update index bounds')
        #also update aa per subplot

        self.param['index'].bounds = (0, len(event.new) - 1)
        self.exposure_str.value = str(self.peptide_measurement.exposure)

        step = 25
        value = int(step*(self.parent.series.cov.end // step + 1))
        self.aa_per_subplot = value# triggers redraw

        #must be uniform
        self.wrap = autowrap(self.parent.series.cov)

        #set index to zero
        self.index = 0

    @param.depends('index', watch=True)
    def _update_index(self):
        self.exposure_str.value = str(self.peptide_measurement.exposure)

    # @property
    # def panel(self):
    #     col = pn.Column(self.param, self.exposure_str)
    #     #p = pn.Param(self.param, widgets={'file': pn.widgets.FileInput}) for exposure
    #     return pn.WidgetBox(col, pn.layout.VSpacer(), css_classes=['widget-box', 'custom-wbox'],
    #                         sizing_mode='stretch_height')


class FittingControl(ControlPanel):
    header = 'Rate Fitting'

    #r_max = param.Number(27, doc='Ceil value for rates', precedence=-1)  # Update this value
    chisq_thd = param.Number(20, doc='Threshold for chi2 to switch to Differential evolution')

    fitting_model = param.Selector(default='Association', objects=['Association', 'Dissociation'])
    do_fit1 = param.Action(lambda self: self._action_fit1())
    block_mode = param.ObjectSelector(default='reduced', objects=['reduced', 'original', 'constant'])

    #todo generate from func signature?
    #block mode reduced params
    max_combine = param.Integer(2, doc='Neighbouring blocks up to and including this size are merged together', precedence=-1)
    max_join = param.Integer(5, doc='Blocks up to and including this size are joined with their smallest neighbour', precedence=-1)

    #constant block params
    block_size = param.Integer(10, doc='Size of the blocks in constant blocks mode', precedence=-1)
    initial_block = param.Integer(5, doc='Size of the initial block in constant block mode', precedence=-1)
    show_blocks = param.Boolean(False, doc='Show how blocks are defined with the current settings', precedence=-1)

    do_fit2 = param.Action(lambda self: self._action_fit2(), constant=True, precedence=-1)

    def __init__(self, parent, **params):
        self.pbar1 = ASyncProgressBar()
        self.pbar2 = ASyncProgressBar()
        super(FittingControl, self).__init__(parent, **params)
                #self.block_column = pn.Column(*[self.param[key] for key in ['max_combine', 'max_join']])
        self.parent.param.watch(self._update_series, ['series'])

    def make_list(self):
        text_f1 = pn.widgets.StaticText(value='Weighted averaging fit (Fit 1)')
        text_f2 = pn.widgets.StaticText(value='Global fit (Fit 2)')

        self._widget_dict.update(text_f1=text_f1, text_f2=text_f2, pbar1=self.pbar1.view, pbar2=self.pbar2.view)
        # parameters = ['r_max', 'fitting_model', 'text_f1', 'chisq_thd', 'do_fit1', 'pbar1', 'text_f2', 'block_mode', 'max_combine',
        #               'max_join', 'show_blocks', 'do_fit2', 'pbar2']
        parameters = ['fitting_model', 'do_fit1', 'pbar1']

        widget_list = list([self._widget_dict[par] for par in parameters])
        return widget_list

    def _update_series(self, *events):
        self.r_max = np.log(1 - 0.98) / -self.parent.series.timepoints[1]  # todo user input 0.98

    async def _fit1_async(self):
        #client = await Client(self.parent.cluster)

        fit_result = await self.parent.fitting.weighted_avg_fit_async(model_type=self.fitting_model.lower(), pbar=self.pbar1)
        print('fit res in async', fit_result)
        self.parent.fit_results['fit1'] = fit_result

        output = fit_result.output
        #todo duplicate code in fit - > method on parent?
        dic = {name: output[name] for name in output.dtype.names}
        dic['y'] = output['rate']  # entry y is by default used for plotting and thresholding
        dic['color'] = np.full_like(output, fill_value=DEFAULT_COLORS['fit1'], dtype='<U7')
        self.parent.sources['fit1'] = ColumnDataSource(dic)

        #trigger plot update
        callback = partial(self.parent.param.trigger, 'sources')
        self.parent.doc.add_next_tick_callback(callback)

        with pn.io.unlocked():
             self.parent.param.trigger('fit_results')  #informs other fittings that initial guesses are now available
             self.pbar1.reset()
             self.param['do_fit1'].constant = False
             self.param['do_fit2'].constant = False

    def _fit1(self):
        fit_result = self.parent.fitting.weighted_avg_fit(model_type=self.fitting_model.lower(), pbar=self.pbar1, chisq_thd=self.chisq_thd)
        self.parent.fit_results['fit1'] = fit_result
        output = fit_result.output
        dic = {name: output[name] for name in output.dtype.names}
        dic['y'] = output['rate']  # entry y is by default used for plotting and thresholding
        dic['color'] = np.full_like(output, fill_value=DEFAULT_COLORS['fit1'], dtype='<U7')
        self.parent.sources['fit1'] = ColumnDataSource(dic)

        self.parent.param.trigger('fit_results')  # Informs TF fitting that now fit1 is available as initial guesses
        self.parent.param.trigger('sources') # Informs listening plots that there is new data available

        self.param['do_fit1'].constant = False
        self.param['do_fit2'].constant = False
        self.pbar1.reset()

    def _action_fit1(self):
        print('fitting 1')
        #todo context manager?
        self.param['do_fit1'].constant = True
        self.param['do_fit2'].constant = True

        if self.parent.cluster:
            print(self.parent.cluster)
            self.parent.doc = pn.state.curdoc
            loop = IOLoop.current()
            loop.add_callback(self._fit1_async)
        else:
            self._fit1()

        # fit_result = self.parent.fitting.weighted_avg_fit(chisq_thd=self.chisq_thd)
        # rates_array = fit_result.get_output(['rate', 'tau', 'tau1', 'tau2', 'r'])

    async def _fit2_async(self):
        fit_result = await self.parent.fitting.blocks_fit_async(self.parent.fit_results['fit1']['rates'], model_type=self.fitting_model.lower(),
                                                                pbar=self.pbar2, **self.fit_kwargs)
        print('fit res in async', fit_result)
        rates_array = fit_result.get_output(['rate', 'tau', 'k1', 'k2', 'r'])
        self.parent.fit_results['fit2'] = {'rates': rates_array, 'fitresult': fit_result}
        callback = partial(self.parent.param.trigger, 'fit_results')
        self.parent.doc.add_next_tick_callback(callback)

        with pn.io.unlocked():
            self.param['do_fit1'].constant = False
            self.param['do_fit2'].constant = False
            self.pbar2.reset()

    def _fit2(self):
        fit_result = self.parent.fitting.blocks_fit(self.parent.fit_results['fit1']['rates'], model_type=self.fitting_model.lower(), **self.fit_kwargs)

        self.parent.fit_results['fit2'] = fit_result
        #todo lines below are repeating code, create func(self, fit_result, name, y=rate) which does this job  (or func on parent controller?)
        output = fit_result.output
        dic = {name: output[name] for name in output.dtype.names}
        dic['y'] = output['rate']
        dic['color'] = np.full_like(output, fill_value=DEFAULT_COLORS['fit2'], dtype='<U7')
        self.parent.sources['fit2'] = ColumnDataSource(dic)
        self.parent.param.trigger('fit_results')  # Informs TF fitting that now fit2 is available as initial guesses
        self.parent.param.trigger('sources')  # Informs listening plots that there is new data available

        #todo all fit buttons should globally constant during fitting (make some kind of global context manager)
        self.param['do_fit1'].constant = False
        self.param['do_fit2'].constant = False
        self.pbar2.reset()

    def _action_fit2(self):
        print('fitting 2')
        #todo context manager
        self.param['do_fit1'].constant = True
        self.param['do_fit2'].constant = True

        if self.parent.cluster:
            self.parent.doc = pn.state.curdoc
            loop = IOLoop.current()
            loop.add_callback(self._fit2_async)
        else:
            self._fit2()

    @property
    def fit_kwargs(self):
        if self.block_mode == 'reduced':
            fit_kwargs = {'block_func': get_reduced_blocks, 'max_combine': self.max_combine, 'max_join': self.max_join}
        elif self.block_mode == 'original':
            fit_kwargs = {'block_func': get_original_blocks}
        elif self.block_mode == 'constant':
            fit_kwargs = {'block_func': get_constant_blocks, 'block_size': self.block_size, 'initial_block': self.initial_block}
        return fit_kwargs

    @property
    def fit_block_edges(self):
        """returns the position of block edges from the current block func"""
        kwargs = self.fit_kwargs
        func = kwargs.pop('block_func')
        all_edges = []
        for series in self.parent.series.split().values():
            blocks = np.array(func(series.cov, **kwargs))
            indices = list(np.cumsum(blocks) - 1)
            edges = [series.cov.start - 0.5] + list(series.cov.r_number[[indices]] + 0.5)
            all_edges += edges

        all_edges = np.array(all_edges)
        return all_edges

    def _clear_block_kwargs(self):
        """removes all block func kwarg widgets from the box"""
        parameters = ['max_combine', 'max_join', 'initial_block', 'block_size']
        for par in parameters:
            try:
                self.box_pop(par)
            except ValueError:
                pass

    @param.depends('block_mode', watch=True)
    def _update_block_mode(self):
        print('block mode updated')
        if self.block_mode == 'reduced':
            self._clear_block_kwargs()
            self.box_insert_after('block_mode', 'max_combine')
            self.box_insert_after('max_combine', 'max_join')
        elif self.block_mode == 'original':
            self._clear_block_kwargs()
        elif self.block_mode == 'constant':
            self._clear_block_kwargs()
            self.box_insert_after('block_mode', 'initial_block')
            self.box_insert_after('initial_block', 'block_size')


class TFFitControl(ControlPanel):
    header = 'TF Single residue fit'

    fitting_type = param.Selector(objects=['Protection Factors', 'Rates'], default='Protection Factors')

    c_term = param.Integer(None, doc='Residue number to which the last amino acid in the sequence corresponds')  # remove
    temperature = param.Number(293.15, doc='Deuterium labelling temperature in Kelvin')
    pH = param.Number(8., doc='Deuterium labelling pH', label='pH')

    stop_loss = param.Number(0.1, bounds=(0, None), doc='Threshold loss difference below which to stop fitting')
    stop_patience = param.Integer(50, bounds=(1, None), doc='Number of epochs where stop loss should be satisfied before stopping')
    learning_rate = param.Number(0.01, bounds=(0, None), doc='Learning rate parameter for optimization')
    epochs = param.Number(10000, bounds=(1, None), doc='Maximum number of epochs (iterations')

    l1_regularizer = param.Number(100, bounds=(0, None), doc='Value for l1 regularizer')
    l2_regularizer = param.Number(0, bounds=(0, None), doc='Value for l2 regularizer')

    do_fit = param.Action(lambda self: self._do_fitting(), constant=True)

    def __init__(self, parent, **params):
        self.pbar1 = ASyncProgressBar()
        super(TFFitControl, self).__init__(parent, **params)

        self.parent.param.watch(self._parent_fit_results_updated, ['fit_results'])
        self.parent.param.watch(self._parent_series_updated, ['series'])

    def make_dict(self):
        kwargs = {name: pn.param.LiteralInputTyped(param.Number(0.)) for name in ['temperature', 'pH', 'stop_loss', 'learning_rate', 'l1_regularizer', 'l2_regularizer']}
        return self.generate_widgets(**kwargs)

    @param.depends('fitting_type', watch=True)
    def _update_fitting_type(self):
        if self.fitting_type == 'Protection Factors':
            self.box_insert_after('fitting_type', 'c_term')
            self.box_insert_after('c_term', 'temperature')
            self.box_insert_after('temperature', 'pH')
        elif self.fitting_type == 'Rates':
            self.box_pop('c_term')
            self.box_pop('temperature')
            self.box_pop('pH')

    def _parent_series_updated(self, *events):
        end = self.parent.series.cov.end
        self.c_term = int(end + 5)

    def _parent_fit_results_updated(self, *events):
        if 'fit1' in self.parent.fit_results:
            self.param['do_fit'].constant = False

    def _do_fitting(self):
        self.param['do_fit'].constant = True

        import pyhdx.fitting_tf as tft
        if self.fitting_type == 'Protection Factors':
            k_int = self.parent.series.cov.calc_kint(self.temperature, self.pH, c_term=self.c_term)
            k_r_number = self.parent.series.cov.sequence_r_number
            k_dict = {'r_number': k_r_number, 'k_int': k_int}

            output_name = 'pfact'
            var_name = 'log_P'
        elif self.fitting_type == 'Rates':
            k_dict = None
            output_name = 'TF_rate'
            var_name = 'log_k'

        initial_result = self.parent.fit_results['fit1'].output

        early_stop = tft.EarlyStopping(monitor='loss', min_delta=self.stop_loss, patience=self.stop_patience)
        result = self.parent.fitting.global_fit(initial_result=initial_result, k_int=k_dict,
                                                learning_rate=self.learning_rate, l1=self.l1_regularizer,
                                                l2=self.l2_regularizer, epochs=self.epochs, callbacks=[early_stop])

        output_dict = {name: result.output[name] for name in result.output.dtype.names}
        output_dict['color'] = np.full_like(result.output, fill_value=DEFAULT_COLORS['pfact'], dtype='<U7')

        output_dict['y'] = 10**output_dict[var_name]
        if self.fitting_type == 'Protection Factors':
            deltaG = constants.R * self.temperature * np.log(output_dict['y'])
            output_dict['deltaG'] = deltaG

        source = ColumnDataSource(output_dict)
        self.parent.sources[output_name] = source
        self.parent.fit_results[output_name] = result

        self.parent.param.trigger('sources')  # dont need to trigger fit_results as its has no relevant watchers
        self.param['do_fit'].constant = False
    #


class FittingQuality(ControlPanel):
    header = 'Fitting Quality'

    peptide_index = param.Number(0, bounds=(0, None))
    x_axis_type = param.Selector(default='Log', objects=['Linear', 'Log'])
    chi_sq = param.Number(0., bounds=(0, None))

    def __init__(self, parent, **param):
        super(FittingQuality, self).__init__(parent, **param)

        self.parent.param.watch(self._series_updated, ['series'])

    def _series_updated(self, *events):

        self.param['peptide_index'].bounds =(0, len(self.parent.series.cov.data))


class ClassificationControl(ControlPanel):
    header = 'Classification'
    num_classes = param.Number(3, bounds=(1, 10), doc='Number of classification classes')
    target = param.Selector(label='Target')
    otsu_thd = param.Action(lambda self: self._action_threshold(), label='Otsu')
    show_thds = param.Boolean(True, label='Show Thresholds')
    values = param.List(precedence=-1)
    colors = param.List(precedence=-1)

    def __init__(self, parent, **param):
        super(ClassificationControl, self).__init__(parent, **param)

        self.values_widgets = []
        for _ in range(self.num_classes - 1):
            self._add_value()

        self.colors_widgets = []
        for _ in range(self.num_classes):
            self._add_color()

        self.param.trigger('values')
        self.param.trigger('colors')
        self.parent.param.watch(self._parent_sources_updated, ['sources'])

    def make_dict(self):
        return self.generate_widgets(num_classes=pn.widgets.Spinner)

    def _parent_sources_updated(self, *events):
        print('sources')
        print("UPDATE")

        objects = list(self.parent.sources.keys())
        self.param['target'].objects = objects

        #set target if its not set already
        if not self.target and objects:
            self.target = objects[-1]

    @param.depends('values', watch=True)
    def _action_threshold(self):
        if self.num_classes > 1 and self.target:
            y_vals = self.parent.sources[self.target].data['y']
            thd_vals = y_vals[~np.isnan(y_vals)]
            thds = threshold_multiotsu(np.log(thd_vals), classes=self.num_classes)
            for thd, widget in zip(thds, self.values_widgets):
                widget.value = np.exp(thd)
        self._do_thresholding()

    def _do_thresholding(self):
        # perhaps we need a class to handle fitting output which has this method
        # yes we do. for all fitting not just fit1
        # alright great. now stop talking to yourself and get back to worK!
        # #quarantine

        # dont do thresholding if the following criteria are met
        if 0 in self.values:
            return
        elif np.any(np.diff(self.values)) < 0:
            return

        y_vals = self.parent.sources[self.target].data['y']
        colors = np.empty(len(y_vals), dtype='U7')

        if self.num_classes == 1:
            colors[:] = self.colors[0]
        else:
            full_thds = [-np.inf] + self.values + [np.inf]
            for lower, upper, color in zip(full_thds[:-1], full_thds[1:], self.colors[::-1]):
                b = (y_vals > lower) & (y_vals <= upper)
                colors[b] = color

       # if 'color' in self.parent.rates.dtype.names:
        print('values', self.values)
        print(self.colors)
        print(colors)
        self.parent.sources[self.target].data['color'] = colors  # this should trigger an update of the graph

    @param.depends('num_classes', watch=True)
    def _update_num_colors(self):
        while len(self.colors_widgets) != self.num_classes:
            if len(self.colors_widgets) > self.num_classes:
                self._remove_color()
            elif len(self.colors_widgets) < self.num_classes:
                self._add_color()
        self.param.trigger('colors')

    @param.depends('num_classes', watch=True)
    def _update_num_values(self):
        while len(self.values_widgets) != self.num_classes - 1:
            if len(self.values_widgets) > self.num_classes - 1:
                self._remove_value()
            elif len(self.values_widgets) < self.num_classes - 1:
                self._add_value()
        print('num classes trigger')
        print(self.values)

        #self._action_threshold()
        self.param.trigger('values')

    def _add_value(self):
        default = 0.0
        self.values.append(default)

        name = 'Threshold {}'.format(len(self.values_widgets) + 1)
        widget = pn.widgets.LiteralInput(name=name, value=default)
        self.values_widgets.append(widget)
        i = len(self.values_widgets) + self.box_index('show_thds')
        self._box.insert(i, widget)
        widget.param.watch(self._value_event, ['value'])

    def _remove_value(self):
        print('remove')
        widget = self.values_widgets.pop(-1)
        self.box_pop(widget)

        self.values.pop()
        print(self.values)

        [widget.param.unwatch(watcher) for watcher in widget.param._watchers]
        del widget

    def _add_color(self):
        try:
            default = DEFAULT_CLASS_COLORS[len(self.colors_widgets)]
        except IndexError:
            default = '#FFFFFF'  #random color?

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

    #todo jslink?
    def _color_event(self, *events):
        print('color event')
        for event in events:
            print(event)
            idx = list(self.colors_widgets).index(event.obj)
            self.colors[idx] = event.new
            c_array = self.parent.sources[self.target].data['color'].copy()
            c_array[c_array == event.old] = event.new
            self.parent.sources[self.target].data['color'] = c_array
        # self.param.trigger('colors')  # i dont think anyone listens to this
        # self.parent.param.trigger('rate_colors')

    def _value_event(self, *events):
        print('value event')
        for event in events:
            idx = list(self.values_widgets).index(event.obj)
            self.values[idx] = event.new
        print(self.values)

        #self._action_threshold()
        self.param.trigger('values')


class FileExportPanel(ControlPanel):
    header = "File Export"
    target = param.Selector(label='Target')

    #todo link this number with the other one
    c_term = param.Integer(0, bounds=(0, None))

    def __init__(self, parent, **param):
        self.export_linear_download = pn.widgets.FileDownload(filename='<no data>', callback=self.linear_export_callback)
        self.pml_script_download = pn.widgets.FileDownload(filename='<no data>', callback=self.pml_export_callback)
        super(FileExportPanel, self).__init__(parent, **param)

        self.parent.param.watch(self._sources_updated, ['sources'])
        self.parent.param.watch(self._series_updated, ['series'])

    def make_list(self):
        self._widget_dict.update(export_linear_download=self.export_linear_download, pml_script_download=self.pml_script_download)
        return super(FileExportPanel, self).make_list()

    def _sources_updated(self, *events):
        objects = list(self.parent.sources.keys())
        self.param['target'].objects = objects

        if not self.target and objects:
            self.target = objects[0]

    def _series_updated(self, *events):
        print('rates updated in fileexportpanel')
        self.c_term = int(self.parent.series.cov.end)
        # #todo centralize this on parent? -> no child controls should hook into main controller
        ## TODO USE .link() function: https://github.com/holoviz/panel/issues/1462

    def _make_pml(self, target):
#        try:
        data_dict = self.parent.sources[target].data
        array = data_dict['y']
        bools = ~np.isnan(array)
        script = colors_to_pymol(data_dict['r_number'][bools], data_dict['color'][bools], c_term=self.c_term)


        return script

#     def pml1_export(self):
#         io = StringIO()
#         script = self._make_pml('fit1')
#         io.write(script)
#         io.seek(0)
#         return io
#
#     def pml2_export(self):
#         io = StringIO()
#         script = self._make_pml('fit2')
#         io.write(script)
#         io.seek(0)
#         return io
#
#     def fit1_export(self):
#         io = StringIO()
# #        print(self.target)
#         print('exporting fit1')
#
#         fit_arr = self.parent.fit_results['fit1']['rates']
#         if 'fit1' in self.parent.rate_colors:
#             colors = self.parent.rate_colors['fit1']
#             export_data = append_fields(fit_arr, 'color', data=colors, usemask=False)
#         else:
#             export_data = fit_arr
#
#         fmt, header = fmt_export(export_data)
#         np.savetxt(io, export_data, fmt=fmt, header=header)
#
#         io.seek(0)
#         return io
#
#     def fit2_export(self):
#         io = StringIO()
#         #        print(self.target)
#         print('exporting fit2')
#
#         fit_arr = self.parent.fit_results['fit2']['rates']
#         if 'fit2' in self.parent.rate_colors:  #todo this should be try/except
#             colors = self.parent.rate_colors['fit2']
#             export_data = append_fields(fit_arr, 'color', data=colors, usemask=False)
#         else:
#             export_data = fit_arr
#
#         fmt, header = fmt_export(export_data)
#         np.savetxt(io, export_data, fmt=fmt, header=header)
#
#         io.seek(0)
#         return io
    @pn.depends('target', watch=True)
    def _update_filename(self):
        self.export_linear_download.filename = self.parent.series.state + '_' + self.target + 'linear.txt'
        self.pml_script_download.filename = self.parent.series.state + '_' + self.target + 'pymol.pml'

    @pn.depends('target')
    def pml_export_callback(self):
        if self.target:
            io = StringIO()
            io.write('# ' + VERSION_STRING + ' \n')
            script = self._make_pml(self.target)
            io.write(script)
            io.seek(0)
            return io
        else:
            return None

    @pn.depends('target')  # param.depends?
    def linear_export_callback(self):
        io = StringIO()
        io.write('# ' + VERSION_STRING + ' \n')

        print(self.target)
        print('exporting')
        if self.target:
            export_dict = {k: v for k, v in self.parent.sources[self.target].data.items() if k != 'y'}
            dtype = [(name, arr.dtype) for name, arr in export_dict.items()]
            export_data = np.empty_like(self.parent.sources[self.target].data['r_number'], dtype=dtype)
            for name, arr in export_dict.items():
                export_data[name] = arr

            fmt, header = fmt_export(export_data)
            np.savetxt(io, export_data, fmt=fmt, header=header)

            io.seek(0)
            return io
        else:
            return None


    # def data_export(self):
    #     io = StringIO()
    #     delimiter = ','
    #
    #     #todo combine these lines into one function?
    #     fmt, header = fmt_export(self.parent.data, delimiter=delimiter, width=0)
    #     np.savetxt(io, self.parent.data, fmt=fmt, header=header, delimiter=delimiter)
    #     io.seek(0)
    #     return io


class OptionsPanel(ControlPanel):
    header = 'Options'

    """panel for various options and settings"""

    #todo this needs to access other panels as well

    link_xrange = param.Boolean(False)

    def __init__(self, parent, **param):
        super(OptionsPanel, self).__init__(parent, **param)
        self.cov_fig_panel = None
        self.rate_fig_panel = None
        self.coverage_ctrl = None

    #
    # def setup(self, fig1, fig2, coverage_ctrl):
    #     self.fig1 = fig1,
    #     self.fig2 = fig2

    @property
    def fig1(self):
        return self.cov_fig_panel.figures[0]

    @property
    def fig2(self):
        return self.rate_fig_panel.figure

    @property
    def enabled(self):
        return self.fig1 is not None and self.fig1 is not None and self.coverage_ctrl is not None

    @param.depends('link_xrange', watch=True)
    def _update_link(self):
        if self.enabled:
            if self.link_xrange:
                self._link()
            else:
                self._unlink()

    def _unlink(self):
        self.coverage_ctrl.param['aa_per_subplot'].constant = False
        self.fig1.x_range.js_property_callbacks.pop('change:start')
        self.fig1.x_range.js_property_callbacks.pop('change:end')

        self.fig2.x_range.js_property_callbacks.pop('change:start')
        self.fig2.x_range.js_property_callbacks.pop('change:end')

    def _link(self):
        step = 25  #todo global config
        value = int(step*(self.parent.series.cov.end // step + 1))
        self.coverage_ctrl.aa_per_subplot = value# triggers redraw
        self.coverage_ctrl.param['aa_per_subplot'].constant = True

        self.fig1.x_range.js_link('start', self.fig2.x_range, 'start')
        self.fig1.x_range.js_link('end', self.fig2.x_range, 'end')

        self.fig2.x_range.js_link('start', self.fig1.x_range, 'start')
        self.fig2.x_range.js_link('end', self.fig1.x_range, 'end')

    # def panel(self):
    #     return pn.WidgetBox(pn.Param(self.param))


class DeveloperPanel(ControlPanel):
    header = 'Developer Options'
    parse = param.Action(lambda self: self._action_load_files())

    def __init__(self, parent, **params):
        self.keys = ['fit1_rates', 'fit1_result', 'fit2_rates', 'fit2_result']
        self.file_selectors = {key: pn.widgets.FileInput() for key in self.keys}
        super(DeveloperPanel, self).__init__(parent, **params)

    def make_list(self):
        return list(self.file_selectors.values()) + [self._widget_dict['parse']]

    def _action_load_files(self):

        for k, fs in self.file_selectors.items():
            if fs.value is not None:
                name = k.split('_')[0]  # fit 1 or fit2
                if 'rates' in k:
                    s_io = StringIO(fs.value.decode('UTF-8'))
                    data = np_from_txt(s_io)
                    self.parent.fit_results[name]['rates'] = data
                elif 'result' in k:
                    b_io = BytesIO(fs.value)
                    result = pickle.load(b_io)
                    self.parent.fit_results[name]['fitresult'] = result
        self.parent.param.trigger('fit_results')

