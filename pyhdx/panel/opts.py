import param
import panel as pn

class Opts(param.Parameterized):

    def __init__(self, **params):
        self._opts = params.pop('opts', {})
        self._excluded_from_opts = list(params.keys())
        super().__init__(**params)

        self.widgets = self.generate_widgets()

    @property
    def panel(self):
        return pn.Column(*self.widgets.values())

    @property
    def opts(self):
        return {**self._opts, **{name: self.param[name] for name in self.param if name not in self._excluded_from_opts}}

    def generate_widgets(self, **kwargs):
        """returns a dict with keys parameter names and values default mapped widgets"""
        #todo base class?

        names = [p for p in self.param if self.param[p].precedence is None or self.param[p].precedence > 1]
        widgets = pn.Param(self.param, show_name=False, show_labels=True, widgets=kwargs)

        return {k: v for k, v in zip(names[1:], widgets)}


class CmapOpts(Opts):

    cmap = param.ObjectSelector(default='jet', objects=['viridis', 'plasma', 'magma', 'jet'], label='Color map')
