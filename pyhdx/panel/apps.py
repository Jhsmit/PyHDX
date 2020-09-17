from pyhdx.panel.template import GoldenElvis, ExtendedGoldenTemplate
from pyhdx.panel.theme import ExtendedGoldenDarkTheme, ExtendedGoldenDefaultTheme
from pyhdx.panel.controller import *
from pyhdx.panel.fig_panels import *
from pyhdx.panel.log import get_default_handler
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


# class DifferenceClassificationControl(ClassificationControl):
#     name = 'ClassificationControl'
#     __name__ = 'ClassificationControl'
#     # log_space = param.Boolean(False, precedence=-1,
#     #                           doc='Boolean to set whether to apply colors in log space or not.')



control_panels = [
    MappingFileInputControl,
    DifferenceControl,
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
    ThdLinearFigure,
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
                             elvis.view(ctrl.figure_panels['ThdLinearFigure']),
                             elvis.view(ctrl.figure_panels['LoggingFigure']),
                         )
                     )
                    )


ctrl.control_panels['ClassificationControl'].log_space = False
ctrl.control_panels['ClassificationControl'].param['log_space'].constant = True


if __name__ == '__main__':
    pn.serve(tmpl, show=True)


