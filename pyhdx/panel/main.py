from pyhdx.panel.controller import Controller
import panel as pn


ctrl = Controller('template', ['asdf'], cluster=None)
ctrl.serve()
