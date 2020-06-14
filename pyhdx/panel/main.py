from pyhdx.panel.controller import Controller
import panel as pn

cluster= '127.0.0.1:59324'
ctrl = Controller('template', ['asdf'], cluster=cluster)
ctrl.serve()
