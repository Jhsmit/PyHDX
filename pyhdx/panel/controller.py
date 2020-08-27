from .log import setup_custom_logger
from .base import ControlPanel, DEFAULT_COLORS, DEFAULT_CLASS_COLORS
from .fig_panels import FigurePanel
from pyhdx.models import PeptideMasterTable, KineticsSeries
from pyhdx.fitting import KineticsFitting
from pyhdx.fileIO import read_dynamx
from pyhdx.support import fmt_export, np_from_txt, \
    autowrap, colors_to_pymol, rgb_to_hex, gen_subclasses

from pyhdx import VERSION_STRING, VERSION_STRING_SHORT


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
from functools import partial
#from .widgets import NumericInput
from bokeh.models import ColumnDataSource, LinearColorMapper, ColorBar
from bokeh.plotting import figure
from collections import namedtuple
#dev only
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
import itertools

from .template import ExtendedGoldenTemplate
from .theme import ExtendedGoldenDarkTheme, ExtendedGoldenDefaultTheme
from .widgets import ColoredStaticText


logger = setup_custom_logger('root')
logger.debug('main message')


HalfLifeFitResult = namedtuple('HalfLifeFitResult', ['output'])


class Controller(param.Parameterized):
    """
    controller for main panels layout
    and has panels for each tabin the main layout

    """

    fit_results = param.Dict({})
    sources = param.Dict({}, doc='Dictionary of ColumnDataSources available for plotting')
    rate_colors = param.Dict({})  #probably not used
    peptides = param.ClassSelector(PeptideMasterTable, doc='Master list of all peptides')
    series = param.ClassSelector(KineticsSeries, doc='Currently selected kinetic series of peptides')
    fitting = param.ClassSelector(KineticsFitting)

    def __init__(self, control_panels, figure_panels, cluster=None, **params):
        super(Controller, self).__init__(**params)
        self.cluster = cluster
        self.doc = pn.state.curdoc

        available_controllers = {cls.__name__: cls for cls in gen_subclasses(ControlPanel)}
        self.control_panels = {name: available_controllers[name](self) for name in control_panels}

        available_figures = {cls.__name__: cls for cls in gen_subclasses(FigurePanel)}
        self.figure_panels = {name: available_figures[name](self) for name in figure_panels}

        #setup options  #todo automate figure out cross dependencies (via parent?)
        self.control_panels['OptionsPanel'].link_xrange = True

    def publish_data(self, name, dic):
        """
        Publish dataset to be available for client figure to plot

        Parameters
        ----------
        name: :obj:`str`
            Name of the dataset
        dic: :obj:`dict`
            Data dictionary, where every key, value pair corresponds to a data column to plot
        """
        source = ColumnDataSource(dic)
        try:  # update existing source
            src = self.sources[name]
            src.data.update(**source.data)
        except KeyError:
            self.sources[name] = source

        self.param.trigger('sources')


class FileInputControl(ControlPanel):
    header = 'Input'

    add_button = param.Action(lambda self: self._action_add(), doc='Add File', label='Add File')
    clear_button = param.Action(lambda self: self._action_clear(), doc='Clear files', label='Clear Files')
    drop_first = param.Integer(1, bounds=(0, None))
    ignore_prolines = param.Boolean(True, constant=True, doc='Set to True to ignore prolines in the sequence')
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
        i = len(self.file_selectors)  # position to insert the new file selector into the widget box
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
    color_map = param.Selector(objects=['jet', 'inferno', 'viridis', 'cividis', 'plasma', 'cubehelix'], default='jet')
    index = param.Integer(0, bounds=(0, 10), doc='Current index of coverage plot in time')

    def __init__(self, parent, **params):
        self.exposure_str = ColoredStaticText(name='Exposure', value='0')  # todo update to some param?

        # We need a reference to color mapper to update it when the cmap changes
        self.color_mapper = LinearColorMapper(palette=self.palette, low=0, high=100)
        self.color_bar = self.get_color_bar()

        super(CoverageControl, self).__init__(parent, **params)
        self.parent.param.watch(self._series_updated, ['series'])

    def make_list(self):
        lst = super(CoverageControl, self).make_list()
        return lst + [self.exposure_str]#, self.color_bar]

    def make_dict(self):
        return self.generate_widgets(index=pn.widgets.IntSlider)

    @param.depends('color_map', watch=True)
    def _update_cbar(self):
        cmap = mpl.cm.get_cmap(self.color_map)
        pal = tuple(mpl.colors.to_hex(cmap(value)) for value in np.linspace(0, 1, 1024, endpoint=True))
        self.color_mapper.palette = pal

    def get_color_bar(self):
        """pn.pane.Bokeh: bokeh pane with empty figure and only a color bar"""
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
        """"""
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
    def colors(self):
        """~class:`np.ndarray`: array of color for each peptide based on their uptake score"""
        cmap = mpl.cm.get_cmap(self.color_map)
        c_rgba = cmap(self.peptide_measurement.data['scores'] / 100)
        c = [mpl.colors.to_hex(color) for color in c_rgba]

        return np.array(c)

    def _series_updated(self, event):
        print('coverage new series update index bounds')

        # must be uniform
        self.wrap = autowrap(self.coverage)
        self.param['index'].bounds = (0, len(event.new) - 1)

        # set index to zero
        self.index = 0
        # self.exposure_str.value = str(self.peptide_measurement.exposure) # this should be triggered

        width = self.coverage.data['end'] - self.coverage.data['start'] # Bars are inclusive, inclusive
        x = self.coverage.data['start'] - 0.5 + (width / 2)
        y = list(itertools.islice(itertools.cycle(range(self.wrap, 0, -1)), len(self.coverage)))
        index = [str(i) for i in range(len(self.coverage))]

        plot_dict = dict(x=x, y=y, width=width, color=self.colors, index=index)
        prop_dict = {name: self.peptide_measurement.data[name] for name in self.peptide_measurement.data.dtype.names}
        dic = {**plot_dict, **prop_dict}

        self.parent.publish_data('coverage', dic)

    @param.depends('wrap', watch=True)
    def _update_wrap(self):
        y = list(itertools.islice(itertools.cycle(range(self.wrap, 0, -1)), len(self.coverage)))
        try:
            self.parent.sources['coverage'].data.update(y=y)
        except KeyError:
            pass

    @param.depends('index', 'color_map', watch=True)
    def _update_colors(self):
        self.exposure_str.value = str(self.peptide_measurement.exposure)  #todo this should be an js_link?
        try:
            tooltip_fields = {field: self.peptide_measurement.data[field] for field in ['scores', 'uptake', 'uptake_corrected']}
            self.parent.sources['coverage'].data.update(color=self.colors, **tooltip_fields)

        except KeyError:
            pass


class FittingControl(ControlPanel):
    header = 'Initial Guesses'
    fitting_model = param.Selector(default='Half-life (λ)', objects=['Half-life (λ)', 'Association', 'Dissociation'])
    do_fit1 = param.Action(lambda self: self._action_fit())

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
        self.parent.publish_data('fit1', dic)

        #trigger plot update
        callback = partial(self.parent.param.trigger, 'sources')
        self.parent.doc.add_next_tick_callback(callback)

        with pn.io.unlocked():
             self.parent.param.trigger('fit_results')  #informs other fittings that initial guesses are now available
             self.pbar1.reset()
             self.param['do_fit1'].constant = False

    def _fit1(self):
        fit_result = self.parent.fitting.weighted_avg_fit(model_type=self.fitting_model.lower(), pbar=self.pbar1, chisq_thd=self.chisq_thd)
        self.parent.fit_results['fit1'] = fit_result
        output = fit_result.output
        dic = {name: output[name] for name in output.dtype.names}
        dic['y'] = output['rate']  # entry y is by default used for plotting and thresholding
        dic['color'] = np.full_like(output, fill_value=DEFAULT_COLORS['fit1'], dtype='<U7')

        self.parent.publish_data('fit1', dic)
        self.parent.param.trigger('fit_results')  # Informs TF fitting that now fit1 is available as initial guesses

        self.param['do_fit1'].constant = False
        self.pbar1.reset()

    def _action_fit(self):
        print('fitting')
        #todo context manager?
        self.param['do_fit1'].constant = True

        if self.fitting_model == 'Half-life (λ)':
            kf = KineticsFitting(self.parent.series)
            output = kf.weighted_avg_t50()
            fit_result = HalfLifeFitResult(output=output)
            dic = {name: output[name] for name in output.dtype.names}
            dic['y'] = output['rate']  # entry y is by default used for plotting and thresholding
            dic['color'] = np.full_like(output, fill_value=DEFAULT_COLORS['half-life'], dtype='<U7')

            self.parent.publish_data('half-life', dic)
            self.parent.fit_results['half-life'] = fit_result

            self.parent.param.trigger('fit_results')  # Informs TF fitting that now fit1 is available as initial guesses

            self.param['do_fit1'].constant = False
        else:

            if self.parent.cluster:
                print(self.parent.cluster)
                self.parent.doc = pn.state.curdoc
                loop = IOLoop.current()
                loop.add_callback(self._fit1_async)
            else:
                self._fit1()


class TFFitControl(ControlPanel):
    header = 'Fitting'

    #fitting_type = param.Selector(objects=['Protection Factors', 'Rates'], default='Protection Factors')
    initial_guess = param.Selector()

    c_term = param.Integer(None, doc='Residue number to which the last amino acid in the sequence corresponds')  # remove
    temperature = param.Number(293.15, doc='Deuterium labelling temperature in Kelvin')
    pH = param.Number(8., doc='Deuterium labelling pH', label='pH')

    stop_loss = param.Number(0.01, bounds=(0, None), doc='Threshold loss difference below which to stop fitting')
    stop_patience = param.Integer(50, bounds=(1, None), doc='Number of epochs where stop loss should be satisfied before stopping')
    learning_rate = param.Number(0.01, bounds=(0, None), doc='Learning rate parameter for optimization')
    epochs = param.Number(100000, bounds=(1, None), doc='Maximum number of epochs (iterations')

    l1_regularizer = param.Number(20, bounds=(0, None), doc='Value for l1 regularizer')
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

    # @param.depends('fitting_type', watch=True)
    # def _update_fitting_type(self):
    #     # (temporarily) removed
    #     if self.fitting_type == 'Protection Factors':
    #         self.box_insert_after('fitting_type', 'c_term')
    #         self.box_insert_after('c_term', 'temperature')
    #         self.box_insert_after('temperature', 'pH')
    #     elif self.fitting_type == 'Rates':
    #         self.box_pop('c_term')
    #         self.box_pop('temperature')
    #         self.box_pop('pH')

    def _parent_series_updated(self, *events):
        end = self.parent.series.cov.end
        self.c_term = int(end + 5)

    def _parent_fit_results_updated(self, *events):
        possible_initial_guesses = ['half-life', 'fit1']
        print(self.parent.fit_results.keys())
        objects = [name for name in possible_initial_guesses if name in self.parent.fit_results.keys()]
        if objects:
            self.param['do_fit'].constant = False

        self.param['initial_guess'].objects = objects
        if not self.initial_guess and objects:
            self.initial_guess = objects[0]

    def _do_fitting(self):
        self.param['do_fit'].constant = True
        import pyhdx.fitting_tf as tft

        kf = KineticsFitting(self.parent.series, temperature=self.temperature, pH=self.pH)
        initial_result = self.parent.fit_results[self.initial_guess].output   #todo initial guesses could be derived from the CDS rather than fit results object
        early_stop = tft.EarlyStopping(monitor='loss', min_delta=self.stop_loss, patience=self.stop_patience)
        result = kf.global_fit_new(initial_result, epochs=self.epochs, learning_rate=self.learning_rate,
                                   l1=self.l1_regularizer, l2=self.l2_regularizer, callbacks=[early_stop])

        output_name = 'pfact'
        var_name = 'log_P'

        output_dict = {name: result.output[name] for name in result.output.dtype.names}
        output_dict['color'] = np.full_like(result.output, fill_value=DEFAULT_COLORS['pfact'], dtype='<U7')
        # output_dict[f'{var_name}_full'] = output_dict[var_name].copy()
        # #todo this should be moved to TFFitresults object (or shoud it?) -> DataObject class (see base) (does coloring)
        # output_dict[var_name][~self.parent.series.tf_cov.has_coverage] = np.nan # set no coverage sections to nan
        output_dict['y'] = 10**output_dict[var_name]
        # if self.fitting_type == 'Protection Factors':
        deltaG = constants.R * self.temperature * np.log(output_dict['y'])
        output_dict['deltaG'] = deltaG

        self.parent.fit_results['fr_' + output_name] = result

        self.parent.publish_data(output_name, output_dict)

        self.param['do_fit'].constant = False
        self.parent.param.trigger('fit_results')
    #


class FittingQuality(ControlPanel):
    header = 'Fitting Quality'

    peptide_index = param.Number(0, bounds=(0, None))
    x_axis_type = param.Selector(default='Log', objects=['Linear', 'Log'])

    def __init__(self, parent, **param):
        super(FittingQuality, self).__init__(parent, **param)

        self.d_uptake = {}  ## Dictionary of arrays (N_p, N_t) with results of fit result model calls
        self.parent.param.watch(self._series_updated, ['series'])
        self.parent.param.watch(self._fit_results_updated, ['fit_results'])

    #@depends(parent.series) will this work? parent should be a classselector param? does that accept subclasses?
    def _series_updated(self, *events):
        self.param['peptide_index'].bounds = (0, len(self.parent.series.cov.data))
        self.d_uptake['uptake_corrected'] = self.parent.series.uptake_corrected.T
        self._update_sources()

    @property
    def fit_timepoints(self):
        time = np.logspace(-2, np.log10(self.parent.series.timepoints.max()), num=250)
        time = np.insert(time, 0, 0.)
        return time

    def _fit_results_updated(self, *events):
        print('fit results updated in fitting quality')
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
            timepoints = self.parent.series.timepoints if name == 'uptake_corrected' else self.fit_timepoints
            dic = {'time': timepoints, 'uptake': array[self.peptide_index, :]}
            self.parent.publish_data(name, dic)

        print('keys', self.parent.sources.keys())


class ClassificationControl(ControlPanel):
    header = 'Classification'
    accepted_sources = ['pfact']  #todo add fit1, half-life etc (or add new data objects)

    target = param.Selector(label='Target')

    mode = param.Selector(default='Discrete', objects=['Discrete', 'Continuous'])#, 'ColorMap'])
    num_colors = param.Number(3, bounds=(1, 10), doc='Number of classification colors')
    #cmap = param.Selector(objects=PRESET_CMAPS)
    otsu_thd = param.Action(lambda self: self._action_otsu(), label='Otsu')
    linear_thd = param.Action(lambda self: self._action_linear(), label='Linear')
    log_space = param.Boolean(True, constant=True, doc='Boolean to set whether to apply colors in log space or not')

    show_thds = param.Boolean(True, label='Show Thresholds')
    values = param.List(precedence=-1)
    colors = param.List(precedence=-1)

    def __init__(self, parent, **param):
        super(ClassificationControl, self).__init__(parent, **param)

        self.values_widgets = []
        for _ in range(self.num_colors - 1):
            self._add_value()

        self.colors_widgets = []
        for _ in range(self.num_colors):
            self._add_color()

        self.param.trigger('values')
        self.param.trigger('colors')
        self.parent.param.watch(self._parent_sources_updated, ['sources'])

    def make_dict(self):
        return self.generate_widgets(num_colors=pn.widgets.Spinner, mode=pn.widgets.RadioButtonGroup)

    def _parent_sources_updated(self, *events):
        print('sources')
        print("UPDATE")

        excluded = ['coverage']
        objects = [key for key in self.parent.sources.keys() if key in self.accepted_sources]
        self.param['target'].objects = list(objects)

        #set target if its not set already
        if not self.target and objects:
            self.target = objects[-1]

        if self.values:
            self._get_colors()

    @property
    def target_array(self):
        """returns the array to calculate colors from"""

        y_vals = self.parent.sources[self.target].data['y']
        return y_vals[~np.isnan(y_vals)]

    def _action_otsu(self):
        if self.num_colors > 1 and self.target:
            #y_vals = self.parent.sources[self.target].data['y']
            #thd_vals = y_vals[~np.isnan(y_vals)]
            func = np.log if self.log_space else lambda x: x
            thds = threshold_multiotsu(func(self.target_array), classes=self.num_colors)
            for thd, widget in zip(thds, self.values_widgets):
                widget.value = np.exp(thd)
        self._get_colors()

    def _action_linear(self):
        if self.log_space:
            thds = np.logspace(np.log(np.min(self.target_array)), np.log(np.max(self.target_array)), num=self.num_colors, endpoint=True, base=np.e)
            for thd, widget in zip(thds, self.values_widgets):
                widget.value = thd

    @param.depends('mode', watch=True)
    def _mode_updated(self):
        if self.mode == 'Discrete':
            self.otsu_thd.constant = False
        elif self.mode == 'Continuous':
            self.otsu_thd.constant = True
        elif self.mode == 'ColorMap':
            self.num_colors = 2
            #todo adjust add/ remove color widgets methods

        self.param.trigger('num_colors')

    @param.depends('values', 'colors', watch=True)
    def _get_colors(self):
        pass

        if 0 in self.values:
            return
        elif np.any(np.diff(self.values)) < 0:
            return

        y_vals = self.parent.sources[self.target].data['y'] # full array with nan entries

        if self.num_colors == 1:
            colors = np.full(len(y_vals), fill_value=self.colors[0], dtype='U7')
            colors[np.isnan(y_vals)] = np.nan
        elif self.mode == 'Discrete':
            full_thds = [-np.inf] + self.values + [np.inf]
            colors = np.full(len(y_vals), fill_value = np.nan, dtype='U7')
            for lower, upper, color in zip(full_thds[:-1], full_thds[1:], self.colors[::-1]):
                b = (y_vals > lower) & (y_vals <= upper)
                colors[b] = color
        elif self.mode == 'Continuous':
            func = np.log if self.log_space else lambda x: x
            vals_sorted = np.sort(func(self.values))
            norm = plt.Normalize(vals_sorted[0], vals_sorted[-1])#, clip=True) currently there is never anythin clipped?
            nodes = norm(vals_sorted)
            cmap = mpl.colors.LinearSegmentedColormap.from_list("custom_cmap", list(zip(nodes, self.colors[::-1])))
            colors_rgba = cmap(norm(func(y_vals)))
            colors = np.array([rgb_to_hex(int(r*255), int(g*255), int(b*255)) for r, g, b, a in colors_rgba])
            colors[np.isnan(y_vals)] = np.nan

        self.parent.sources[self.target].data['color'] = colors  # this triggers an update of the graph

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
        print('num classes trigger')
        print(self.values)

        #self._action_otsu()
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

    #todo jslink? (no also link to coloring function)
    def _color_event(self, *events):

        for event in events:
            idx = list(self.colors_widgets).index(event.obj)
            self.colors[idx] = event.new

        self.param.trigger('colors')
        print('color event')
        # for event in events:
        #     print(event)
        #     idx = list(self.colors_widgets).index(event.obj)
        #     self.colors[idx] = event.new
        #     c_array = self.parent.sources[self.target].data['color'].copy()
        #     c_array[c_array == event.old] = event.new
        #     self.parent.sources[self.target].data['color'] = c_array
        # self.param.trigger('colors')  # i dont think anyone listens to this
        # self.parent.param.trigger('rate_colors')

    def _value_event(self, *events):
        print('value event')
        for event in events:
            idx = list(self.values_widgets).index(event.obj)
            self.values[idx] = event.new
        print(self.values)

        #self._action_otsu()
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
    @property
    def export_dict(self):
        return self.parent.sources[self.target].data

    @pn.depends('target', watch=True)
    def _update_filename(self):
        self.export_linear_download.filename = self.parent.series.state + '_' + self.target + '_linear.txt'
        if 'r_number' in self.export_dict.keys():
            self.pml_script_download.filename = self.parent.series.state + '_' + self.target + '_pymol.pml'
            self.pml_script_download.disabled = False
            print(self.pml_script_download.disabled)

        else:
            self.pml_script_download.filename = 'Not Available'
            self.pml_script_download.disabled = True

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
            export_dict = {k: np.array(v) for k, v in self.parent.sources[self.target].data.items() if k != 'y'}  #todo generalize export
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

    link_xrange = param.Boolean(False)

    def __init__(self, parent, **param):
        super(OptionsPanel, self).__init__(parent, **param)

    @property
    def enabled(self):
        return self.master_figure is not None and self.client_figures is not None

    @param.depends('link_xrange', watch=True)
    def _update_link(self):
        print('link_xrangek')
        print(self.enabled)
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


class ProteinViewControl(ControlPanel):
    header = 'Protein Viewer'
    accepted_sources = ['pfact']  #todo add fit1, half-life etc (or add new data objects) (refactor to dataset/object)
    #also linked sources dutn work yet
    target_dataset = param.Selector()

    input_option = param.Selector(default='Upload File', objects=['Upload File', 'RCSB PDB'])
    rcsb_id = param.String()
    load_structure = param.Action(lambda self: self._load_structure())
    no_coverage = param.Color(default='#8c8c8c')
    representation = param.Selector(default='cartoon',
                                    objects=['ball+stick', 'backbone', 'ball+stick', 'cartoon', 'hyperball', 'licorice',
                                             'ribbon', 'rope', 'spacefill', 'surface'])
    spin = param.Boolean(default=False)

    def __init__(self, parent, **params):
        self.file_input = pn.widgets.FileInput(accept='.pdb')
        super(ProteinViewControl, self).__init__(parent, **params)

        self.parent.param.watch(self._parent_sources_updated, ['sources'])
        self.input_option = 'RCSB PDB'

    def make_list(self):
        lst = super().make_list()
        lst.pop(2)  # Remove RCSB ID input field
        lst.insert(2, self.file_input)  # add File input widget
        return lst

    def _parent_sources_updated(self, *events):
        objects = [key for key in self.parent.sources.keys() if key in self.accepted_sources]
        self.param['target_dataset'].objects = list(objects)

    @param.depends('input_option', watch=True)
    def _update_input_option(self):
        if self.input_option == 'Upload File':
            self.box_pop('rcsb_id')
            self.box_insert_after('input_option', self.file_input)
        elif self.input_option == 'RCSB PDB':
            self.box_pop(self.file_input)
            self.box_insert_after('input_option', 'rcsb_id')

    # def _load_structure(self):
    #     if self.input_option == 'Upload File':
    #         if self.file_widget.value:
    #             string = self.file_widget.value.decode()
    #             self.ngl_html.pdb_string = string
    #         else:
    #             pass

        elif self.input_option == 'RCSB PDB':
            self.ngl_html.rcsb_id = self.rcsb_id