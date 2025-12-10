from pyhdx.web.controllers import PeptideFileInputControl
from pyhdx.web.main_controllers import MainController
import panel as pn

main = MainController(control_panels=[(PeptideFileInputControl, {})])


input_ctrl = main.control_panels["PeptideFileInputControl"]

# pn.serve(input_ctrl.panel)


input_ctrl.panel.servable()
