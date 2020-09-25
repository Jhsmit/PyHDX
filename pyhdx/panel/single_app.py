from pyhdx.panel.template import GoldenElvis, ExtendedGoldenTemplate
from pyhdx.panel.theme import ExtendedGoldenDarkTheme, ExtendedGoldenDefaultTheme
from pyhdx.panel.controller import *
from pyhdx.panel.fig_panels import *
from pyhdx.panel.log import get_default_handler
from pyhdx.panel.apps import SingleValueFigure

import sys
from pyhdx import VERSION_STRING_SHORT
import panel as pn

DEBUG = False


class DifferenceFileExportControl(FileExportControl):
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

        r_max = int(np.nanmax(self.export_dict['r_number'])) + 5
        if self.c_term < r_max:
            self.c_term = r_max


class SingleControl(ControlPanel):
    header = 'Datasets'

    dataset = param.Selector(doc='ds1')
    dataset_name = param.String()
    quantity = param.Selector(doc="Select a quantity to plot (column from input txt file)")

    add_dataset = param.Action(lambda self: self._action_add_dataset(),
                                  doc='Click to add this comparison to available comparisons')
    dataset_list = param.ListSelector(doc='Lists available comparisons')
    remove_dataset = param.Action(lambda self: self._action_remove_comparison())

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


class SingleFigure(BinaryComparisonFigure):
    title = 'Values'
    y_label = 'Value'


control_panels = [
    MappingFileInputControl,
    SingleControl,
    ClassificationControl,
    # 'FileExportControl',
    ProteinViewControl,
    DifferenceFileExportControl,
    OptionsControl,
    DeveloperControl
]

if DEBUG:
    control_panels.append('DeveloperControl')

figure_panels = [
    #'LogLinearFigure',
    BinaryComparisonFigure,
    ProteinFigure,
    LoggingFigure
]

elvis = GoldenElvis(ExtendedGoldenTemplate, ExtendedGoldenDarkTheme, title=VERSION_STRING_SHORT)
cluster = '127.0.0.1:52123'
ctrl = ComparisonController(control_panels, figure_panels, cluster=cluster)
ctrl.logger.addHandler(get_default_handler(sys.stdout))
tmpl = elvis.compose(ctrl.control_panels.values(),
                     elvis.column(
                         elvis.stack(
                             elvis.view(ctrl.figure_panels['ProteinFigure'])
                         ),
                         elvis.row(
                             elvis.stack(
                                elvis.view(ctrl.figure_panels['BinaryComparisonFigure']),
                             ),
                             elvis.view(ctrl.figure_panels['LoggingFigure']),
                         )
                     )
                    )


ctrl.control_panels['ClassificationControl'].log_space = False
#ctrl.control_panels['ClassificationControl'].param['log_space'].constant = True

# from pyhdx.support import np_from_txt
# import os
# directory = r'C:\Users\jhsmi\pp\pyHDX_paper\v3\secb_comparison\fit'
# array1 = np_from_txt(os.path.join(directory, 'ecSecB_10_fitted_tf_pfact.txt'))
# array2 = np_from_txt(os.path.join(directory, 'ecSecB_20_fitted_tf_pfact.txt'))
# array3 = array2.copy()
# array3['log_P'] += np.random.rand()
#
#
# ctrl.datasets['ds_10'] = array1
# ctrl.datasets['ds_20'] = array2
# ctrl.datasets['ds_30'] = array3
#
# ctrl.param.trigger('datasets')
#
# ctrl.control_panels['ClassificationControl'].log_space = False
# ctrl.control_panels['ProteinViewControl'].rcsb_id = '1qyn'


if __name__ == '__main__':
    pn.serve(tmpl, show=False, websocket_origin=['pyhdx.loca.lt', 'localhost:51337', '9ce0a771395c.ngrok.io'], port=51337)


