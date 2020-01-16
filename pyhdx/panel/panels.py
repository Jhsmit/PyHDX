import os
import panel as pn
import numpy as np
from bokeh.models.widgets import Button as BKButton
from bokeh.models import CustomJS, ColumnDataSource
from io import StringIO
from pyhdx import PeptideCSVFile


class HDXPanel(object):
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



