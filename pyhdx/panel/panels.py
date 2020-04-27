import matplotlib
matplotlib.use('agg')

import os
import panel as pn
import numpy as np
from bokeh.models.widgets import Button as BKButton
from bokeh.models import CustomJS, ColumnDataSource, LabelSet
from bokeh.plotting import figure
from bokeh.layouts import column, row
from io import StringIO
from pyhdx import PeptideCSVFile, KineticsSeries
from pyhdx.fitting import fit_kinetics
from pyhdx.plot import make_kinetics_figure, _bokeh_coverage
from pyhdx.support import get_reduced_blocks, get_constant_blocks
from pyhdx.fileIO import read_dynamx
import param
from collections import namedtuple
from numpy.lib.recfunctions import stack_arrays

from matplotlib.collections import LineCollection
import matplotlib as mpl
import matplotlib.pyplot as plt


import logging
logger = logging.getLogger('pyhdx')


class PanelBase(param.Parameterized):
    """base class for mixin panels"""

    position = ''

    @property
    def panel(self):
        return None


class FileInputPanel(PanelBase):
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
        super(FileInputPanel, self).__init__(**params)
        self.parent = parent
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
    def control_panel(self):
        params = pn.panel(self.param)
        col = pn.Column(*[self.file_selectors_column, *params[1:]])

        return pn.WidgetBox(col, pn.layout.VSpacer(), css_classes=['widget-box', 'custom-wbox'], sizing_mode='stretch_height')


class RateConstantPanel(PanelBase):
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
        super(RateConstantPanel, self).__init__(**params)
        raise DeprecationWarning('this is all old crap')
        self.parent = parent

        self.figure = figure(y_axis_type="log")
        self.bk_pane = pn.pane.Bokeh(self.figure, sizing_mode='stretch_both')

        self.block_column = pn.Column(*[self.param[key] for key in ['max_combine', 'max_join']])

        self.parent.param.watch(self._renew, ['rates'])
        self.parent.param.watch(self._update, ['series'])

    def _update(self, event):
        #rerender plot because of new series

        source = ColumnDataSource({name: self.parent.rates[name] for name in self.parent.rates.dtype.names})

        self.figure.circle(x='r_number', y='fit1', legend_label='Fit 1', source=source)
        self.figure.circle(x='r_number', y='fit2', legend_label='Fit 2', source=source, color='red')

        #self.figure.circle(x='r_number', y='fit1_r1', legend_label='Fit 1 r1', source=source, color='green')
        #self.figure.circle(x='r_number', y='fit1_r2', legend_label='Fit 1 r2', source=source, color='yellow')

        self.bk_pane.param.trigger('object')

    def _renew(self, event):
        print('rates array update, renew')

        self.r_max = np.log(1 - 0.98) / - self.parent.series.times[1]  # KEEP THIS

        new_dict = {name: self.parent.rates[name] for name in self.parent.rates.dtype.names}
        for renderer in self.figure.renderers:
            renderer.data_source.data.update(new_dict)

        self.bk_pane.param.trigger('object')

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
        self._renew(None)  #manual trigger

        self.param['do_fit1'].constant = False
        self.param['do_fit2'].constant = False

    def _action_fit2(self):
        print('fitting 2')
        #todo context manager
        self.param['do_fit1'].constant = True
        self.param['do_fit2'].constant = True

        r_number, fit2_rate = self.parent.fitting.lsq_fit_blocks(self.parent.rates, **self.fit_kwargs)
        self.parent.rates['fit2'] = fit2_rate
        self._renew(None)  # manual trigger

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
    def view_panel(self):
        return self.bk_pane

    @property
    def control_panel(self):
        par1 = ['chisq_thd', 'r_max', 'do_fit1', 'block_mode']
        par2 = ['do_fit2']

        pars = [self.param[key] for key in par1] + [self.block_column] + [self.param[key] for key in par2]
        return pn.WidgetBox(*pars)


class CoveragePanel(PanelBase):
    wrap = param.Integer(25, bounds=(0, None), doc='Number of peptides vertically before moving to the next row') # todo auto?
    aa_per_subplot = param.Integer(100, label='Amino acids per subplot')
    #color = param.Boolean(False, label='Color')

    update = param.Action(lambda self: self._update(), label='Update')

    labels = param.Boolean(False, label='Labels')

    index = param.Integer(0, bounds=(0, None), doc='Current index of coverage plot in time')


    # next_btn = param.Action(lambda self: self._action_next(), label='Next')
    # prev_btn = param.Action(lambda self: self._action_prev(), label='Previous')

    def __init__(self, parent, **params):
        super(CoveragePanel, self).__init__(**params)
        raise DeprecationWarning('will be removed')
        self.parent = parent

        self.figures = [figure()]
        self.layout = column(*self.figures, sizing_mode='stretch_both')
        self.label_set = LabelSet()
        self.bk_pane = pn.pane.Bokeh(self.layout, sizing_mode='stretch_both')

        #self.exposure_str = pn.pane.Str('Exposure: ')
        self.exposure_str = pn.widgets.StaticText(name='Exposure', value='A string')
        self.parent.param.watch(self._update, ['series'])


    def _render_figure(self):
        layout, figs, label_set = _bokeh_coverage(self.peptide_measurement, self.wrap, self.aa_per_subplot)
        label_set.visible = self.labels
        return layout, figs, label_set

    def _update(self, *events):
        self.param['index'].bounds = (0, len(self.parent.series) - 1)  #THIS HAS TO STAY
        self.exposure_str.value = str(self.peptide_measurement.exposure)  # THIS A SWELL

        self.layout, self.figures, self.label_set = self._render_figure()
        self.bk_pane.object = self.layout

        print('series change')

    @property
    def peptide_measurement(self):
        return self.parent.series[self.index]

    def _get_color(self):
        #todo move function to pyhdx.plot
        #todo also update pallette in that function

        cmap = mpl.cm.get_cmap('jet')
        c_rgba = cmap(self.peptide_measurement.data['scores'] / 100)
        c = [mpl.colors.to_hex(color) for color in c_rgba]

        return list(c)

    @param.depends('index', watch=True)
    def _update_index(self):
        print('index updated', self.index)
        color = self._get_color()
        print(self.peptide_measurement.exposure)

        self.exposure_str.value = str(self.peptide_measurement.exposure)  # THIS HAS TO STAY
        for fig in self.figures:
            fig.renderers[0].data_source.data.update({'c': color})
        self.bk_pane.param.trigger('object')

    @param.depends('labels', watch=True)
    def _update_labels(self):
        self.label_set.visible = self.labels
        self.bk_pane.param.trigger('object')

    @property
    def view_panel(self):
        return self.bk_pane

    @property
    def control_panel(self):
        #todo this is standard, move to base class
        #params = pn.panel(self.param)
        col = pn.Column(self.param, self.exposure_str)
        #p = pn.Param(self.param, widgets={'file': pn.widgets.FileInput}) for exposure
        return pn.WidgetBox(col, pn.layout.VSpacer(), css_classes=['widget-box', 'custom-wbox'],
                            sizing_mode='stretch_height')


class ClassificationPanel(PanelBase):
    num_classes = param.Integer(3, doc='Number of rate constant classes')
    otsu_thd = param.Action(lambda self: self._action_otsu)

    def __init__(self, **params):
        super(ClassificationPanel, self).__init__(**params)




class HDXBase(param.Parameterized):
    file = param.FileSelector()
    drop_first = param.Integer(default=1, bounds=(0, None))
    load_button = param.Action(lambda self: self._action_load(), doc='Load')

    control_state = param.Selector(doc='State for the control condition', label='Control State')
    control_exposure = param.Selector(doc='Exposure for control condition', label='Control exposure')

    parse_button = param.Action(lambda self: self._action_parse(), doc="Parse")

    exp_state = param.Selector(doc='State for selected experiment', label='Experiment State')
    exp_times = param.ListSelector(default=[], objects=[''])

    def __init__(self, **params):
        super(HDXBase, self).__init__(**params)
        raise DeprecationWarning('will also be removed')
        self.peptide_file = None  # PeptideCSVFile object
        self.pm_dict = {}  # Dictionary of PeptideMeasurements object

    def _action_load(self):
        s = StringIO(self.file.decode('UTF-8'))
        self.peptide_file = PeptideCSVFile(s, drop_first=self.drop_first)

        states = list(np.unique(self.peptide_file.data['state']))
        self.param['control_state'].objects = states
        self.control_state = states[0]

    def _action_parse(self):
        self.pm_dict = self.peptide_file.return_by_name(self.control_state, self.control_exposure)

        states = list(np.unique([v.state for v in self.pm_dict.values()]))
        self.param['exp_state'].objects = states
        self.exp_state = states[0]

    @param.depends('exp_state', watch=True)
    def _update_exp_exposure(self):
        exposures = list([v.exposure for v in self.pm_dict.values() if v.state == self.exp_state])
        self.param['exp_times'].objects = exposures
        self.exp_times = exposures

    @param.depends('control_state', watch=True)
    def _update_control_exposure(self):
        b = self.peptide_file.data['state'] == self.control_state
        data = self.peptide_file.data[b]
        exposures = list(np.unique(data['exposure']))
        self.param['control_exposure'].objects = exposures
        self.control_exposure = exposures[0]

    @param.output(pm_dict=param.Dict)
    def output(self):
        d_states = {k: v for k, v in self.pm_dict.items() if v.state == self.exp_state and v.exposure in self.exp_times}
        return d_states

    def panel(self):
        p = pn.Param(self.param, widgets={'file': pn.widgets.FileInput})
        return p


EmptyResult = namedtuple('EmptyResult', ['chi_squared', 'params'])


class HDXKinetics(param.Parameterized):
    pm_dict = param.Dict(precedence=-1)
    chi_squared_max = param.Number(default=20, bounds=(0, None), label='Maximum chi squared',
                                   doc='Maximum value of chi squared below which DifferentialEvolution is used')
    rate_max = param.Number(default=100, bounds=(0, None), label='Maximum rate',
                            doc='Maximum rate for fitted rate constants (1/min')
    fitting_button = param.Action(lambda self: self._action_fitting(), doc='Fit', label='Do Fitting')
    fitting_progress = param.Number(default=0, bounds=(0, 100))

    modeling_type = param.ObjectSelector(default='Otsu', objects=['Otsu', 'HMM'])
    num_classes = param.Integer(default=3, bounds=(2, 10))

    update = param.Action()

    def __init__(self, **params):
        raise DeprecationWarning()

        super(HDXKinetics, self).__init__(**params)
        s = self.pm_dict[next(iter(self.pm_dict))]  # First element in dictionary
        for v in self.pm_dict.values():
            assert np.all(s.cs == v.cs), 'Not all entries in the selected data series are equal'
            assert np.all(s.state == v.state), 'Not all entries in the selected data series are equal'
            assert np.all(
                s.data['sequence'] == v.data['sequence']), 'Not all entries in the selected data series are equal'

        sorted_dict = {k: v for k, v in sorted(self.pm_dict.items(), key=lambda item: item[1].exposure)}
        self.times = np.array([v.exposure for v in sorted_dict.values()])  #units minutes

        # Array of weighted avarage scores with rows equal to numer of time points
        scores_2d = np.stack([v.scores_average for v in sorted_dict.values()])

        # Normalized to 100 array of scores
        self.scores_norm = 100 * (scores_2d / scores_2d[-1, :][np.newaxis, :])

        # Cumalitive sum of residue number for blocks of unique residues
        self.cumsum = s.cs

        # Number of residues per unique block (np.cumsum to get cumsum)
        self.counts = s.counts

        # residue number of start of first peptide
        self.start = s.start

        # residue index
        self.r_number = np.arange(s.start, s.stop + 1)

        er = EmptyResult(np.nan, {k: np.nan for k in ['tau1', 'tau2', 'r']})
        self.results = [er for i in self.r_number]  # List of final results param?

        # move plotting logic to mixin / module?
        self.fig, self.axes = self.get_figure()
        self.mpl = pn.pane.Matplotlib(self.fig, dpi=100)

        # reimplement with download widget coming in panel 0.8.0
        self.save_btn = BKButton(label='Save', button_type='success')
        self.source = ColumnDataSource(data=dict())
        self.save_btn.js_on_click(
            CustomJS(args=dict(source=self.source),
                     code=open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "download.js")).read()))

        self.hdx_viewer_btn = BKButton(label='HDX Viewer export', button_type='success')
        self.hdx_source = ColumnDataSource(data=dict())
        self.hdx_viewer_btn.js_on_click(CustomJS(args=dict(source=self.hdx_source, filename='HDX_export.csv'),
                        code=open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "txt_download.js")).read()))

        self.rate_max = - (np.log(1 - 0.95) / self.times[1])  # assuming first timepoint is zero (try/except?)

    def _action_fitting(self):
        # print('fitting')

        self.param['fitting_button'].constant = True
        i = 0
        results = []
        # todo this needs to be threaded (dummy parameters)?

        for n, j in enumerate(self.cumsum):
            arr = self.scores_norm[:, i:j]
            i = j
            print(n)
            if np.all(np.isnan(arr)):
                result = EmptyResult(np.nan, {k: np.nan for k in ['tau1', 'tau2', 'r']})
            else:
                assert np.all(np.std(arr, axis=1) < 1e-10)
                d = arr[:, 0]
                result = fit_kinetics(self.times, d, self.chi_squared_max)

            progress = int(100 * (n / len(self.cumsum)))
            self.fitting_progress = progress
            results.append(result)

        self.results = np.repeat(results, self.counts)
        self.fitting_progress = 100

        self._update_hdx_viewer_export()

        self.param['fitting_button'].constant = False
        self.param.trigger('update')

        # self.param['download_button'].constant = False

    @param.depends('update', watch=True)
    def _update_export_data(self):
        export_data = self.export_data
        data_dict = {name: export_data[name] for name in export_data.dtype.names}
        self.source.data = data_dict

    @param.depends('update', watch=True)
    def _update_hdx_viewer_export(self):
        s = 'Residues,k,Ln(k)\n'
        rates = np.clip(self.rate, None, self.rate_max)
        for r, rate in zip(self.r_number, rates):
            if np.isnan(rate):
                continue
            line = '{},{:.3},{:.3}'.format(r, rate, np.log10(rate))
            s += line
            s += '\n'

        self.hdx_source.data = {'text': [s]}

    @property
    def export_data(self):
        """returns data for saving/export"""
        dtype = [('position', int), ('tau', float), ('rate', float), ('chi_squared', float)]
        data = np.empty(len(self.r_number), dtype=dtype)
        data['position'] = self.r_number

        rate = np.clip(self.rate, None, self.rate_max)
        data['tau'] = 1/rate
        data['rate'] = rate
        data['chi_squared'] = self.chi_squared

        return data

    @property
    def tau(self):
        return np.array(
            [res.params['r'] * res.params['tau1'] + (1 - res.params['r']) * res.params['tau2'] for res in self.results])

    @property
    def rate(self):
        return 1 / self.tau

    @property
    def chi_squared(self):
        return np.array([res.chi_squared for res in self.results])

    # @property
    # def r_number(self):
    #     """residue number array"""
    #     # TODO check +1 or not
    #     return np.arange(len(self.tau)) + self.start

    def panel(self):
        par = pn.Param(self.param, widgets={
            'fitting_progress': {'type': pn.widgets.Progress, 'sizing_mode': 'stretch_both'}})

        col = pn.Column(par, self.save_btn, self.hdx_viewer_btn)
        row = pn.Row(col, self.mpl)

        return row

    @param.depends('modeling_type', 'num_classes', watch=True)
    def _do_modeling(self):
        if self.modeling_type == 'Otsu':
            pass
        elif self.modeling_type == 'HMM':
            pass

    @property
    def _lc_data(self):
        """data for linecollection at max rate"""
        return [np.column_stack([self.r_number, self.rate_max*np.ones_like(self.r_number)])]

    def get_figure(self):
        """returns matplotlib figure for visualization of kinetics"""

        fig, (ax1, ax2, cbar_ax) = make_kinetics_figure(self.pm_dict)

        lc = LineCollection(self._lc_data, colors='r')
        ax2.add_collection(lc)

        return fig, (ax1, ax2, cbar_ax)

    @param.depends('rate_max', watch=True)
    def update_rate_max(self):
        lc = self.axes[1].collections[0]
        lc.set_segments(self._lc_data)
        self.update_figure()

    @param.depends('update', watch=True)
    def update_figure(self):
        """add rates to figure"""
        ax = self.axes[1]
        line = ax.lines[0]
        line.set_ydata(self.export_data['rate'])
        #ax.update_datalim(line.get_datalim(ax.transData))
        ax.relim()


        ax.autoscale_view()
        #ax.plot(self.r_number, self.export_data['rate'], marker='.', linestyle='', color='k')
        self.mpl.param.trigger('object')



class HDXBaseDep(param.Parameterized):
    file = param.FileSelector()
    drop_first = param.Integer(default=1, bounds=(0, None))
    load_button = param.Action(lambda self: self._action_load(), doc='Load')

    control_state = param.Selector(doc='State for the control condition', label='Control State')
    control_exposure = param.Selector(doc='Exposure for control condition', label='Control exposure')

    parse_button = param.Action(lambda self: self._action_parse(), doc="Parse")

    exp_state = param.Selector(doc='State for selected experiment', label='Experiment State')
    exp_times = param.ListSelector(default=[], objects=[''])

    def __init__(self, **params):
        super(HDXBaseDep, self).__init__(**params)
        raise DeprecationWarning()

        self.peptide_file = None
        self.pm_dict = {}  # Dictionary of PeptideMeasurements object

    def _action_load(self):
        s = StringIO(self.file.decode('UTF-8'))
        self.peptide_file = PeptideCSVFile(s, drop_first=self.drop_first)

        states = list(np.unique(self.peptide_file.data['state']))
        self.param['control_state'].objects = states
        self.control_state = states[0]

    def _action_parse(self):
        self.pm_dict = self.peptide_file.return_by_name(self.control_state, self.control_exposure)
        states = list(np.unique([v.state for v in self.pm_dict.values()]))
        self.param['exp_state'].objects = states
        self.exp_state = states[0]

    @param.depends('exp_state', watch=True)
    def _update_exp_exposure(self):
        exposures = list([v.exposure for v in self.pm_dict.values() if v.state == self.exp_state])
        self.param['exp_times'].objects = exposures
        self.exp_times = exposures

    @param.depends('control_state', watch=True)
    def _update_control_exposure(self):
        b = self.peptide_file.data['state'] == self.control_state
        data = self.peptide_file.data[b]
        exposures = list(np.unique(data['exposure']))
        self.param['control_exposure'].objects = exposures
        self.control_exposure = exposures[0]

    def panel(self):
        p = pn.Param(self.param, widgets={'file': pn.widgets.FileInput},
                     expand_layout=pn.Row, expand_button=True, expand=True)
        return p


class HDXPanel(param.Parameterized):
    """deprecated"""


    control_state = param.Selector(doc='State for the control condition')
    control_exposure = param.Number(bounds=(0, None), doc='Exposure for control condition')
    exp_state = param.Selector(doc='State for selected experiment')
    exp_times = param.List()

    def __init__(self):
        raise DeprecationWarning()
        self.pm_dict = {}
        self.source = ColumnDataSource(data=dict())

        self.file_input = pn.widgets.FileInput()
        self.button = pn.widgets.Button(name='Go!', button_type='primary')
        self.text = pn.widgets.TextInput(value='Ready')

        self.button.on_click(self.load_file)
        self.i = 0

        self.select_state = pn.widgets.Select(name='State')
        self.select_exposure = pn.widgets.Select(name='Exposure')
        self.parse_button = pn.widgets.Button(name='Parse')
        self.parse_button.on_click(self.parse)

        self.multi_select = pn.widgets.MultiSelect(name='Datasets', height=250)

        self.process_button = pn.widgets.Button(name='Process', button_type='default')
        self.process_button.on_click(self.process_dld)
        self.save_btn = BKButton(label='Save', button_type='success')
        self.save_btn.js_on_click(
            CustomJS(args=dict(source=self.source), code=open(os.path.join(os.path.abspath(''), "download.js")).read()))

    def process(self, event):
        home = os.path.expanduser('~')

        for k, d in self.pm_dict.items():
            x = np.arange(d.stop + 1) + 1
            y = np.empty_like(x, dtype=float)
            y.fill(np.nan)
            y[d.start - 1:d.stop] = d.scores_average
            out = np.column_stack((x, y))
            np.savetxt(os.path.join(home, k + '.txt'), out, fmt=['%i', '%f'])

    def process_dld(self, event):
        """
        Process and place resulting output in self.source for downloading
        :param event:
        :return:
        """
        # todo fix multi select selection
        selected = self.multi_select.value
        dtype = [(k, float) for k in selected] + [('position', int)]
        size = list(self.pm_dict.values())[0].stop + 1
        out = np.empty(size, dtype=dtype)
        out['position'] = np.arange(size) + 1
        for k in selected:
            d = self.pm_dict[k]
            y = np.empty(size, dtype=float)
            y.fill(np.nan)
            y[d.start - 1:d.stop] = d.scores_average
            out[k] = y

        self.out = out
        data = {name: out[name] for name in out.dtype.names}
        self.source.data = data

    # CustomJS(args=dict(source=cds), code=open(os.path.join(os.path.abspath(''), "download.js")).read())

    def parse(self, event):
        self.pm_dict = self.pf.return_by_name(self.select_state.value, self.select_exposure.value)

        values = list(self.pm_dict.keys())
        self.multi_select.options = values
        self.multi_select.value = values

    def load_file(self, event):
        self.text.value = 'Clicked {0} times'.format(self.button.clicks)
        self.i += 1
        print('load file')
        s = StringIO(self.file_input.value.decode('UTF-8'))
        self.pf = PeptideCSVFile(s)

        states = list(np.unique(self.pf.data['state']))
        self.select_state.options = states

        exposures = list(np.unique(self.pf.data['exposure']))
        self.select_exposure.options = exposures

    @property
    def panel(self):
        row1 = pn.Row(self.file_input, self.button, self.text)
        row2 = pn.Row(self.select_state, self.select_exposure, self.parse_button)
        row3 = pn.Row(pn.Spacer(sizing_mode='stretch_both'), self.process_button, self.save_btn,
                      pn.Spacer(sizing_mode='stretch_both'), )

        panel = pn.Column(row1, row2, self.multi_select, row3)

        return panel



