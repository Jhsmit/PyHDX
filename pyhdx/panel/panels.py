import os
import panel as pn
import numpy as np
from bokeh.models.widgets import Button as BKButton
from bokeh.models import CustomJS, ColumnDataSource

from io import StringIO
from pyhdx import PeptideCSVFile
from pyhdx.fitting import fit_kinetics
from pyhdx.plot import make_kinetics_figure
import param
from collections import namedtuple

from matplotlib.collections import LineCollection
#matplotlib.use('agg')


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
        self.mpl = pn.pane.Matplotlib(self.fig)

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
        # todo this needs to be threaded

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



