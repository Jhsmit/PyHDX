from pyhdx.panel.controller import Controller
import panel as pn

cluster= '127.0.0.1:52123'
ctrl = Controller('template', ['asdf'], cluster=cluster)
ctrl.serve()
