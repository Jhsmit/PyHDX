from pyhdx.web.apps import main_app

ctrl, tmpl = main_app()

tmpl.servable()
