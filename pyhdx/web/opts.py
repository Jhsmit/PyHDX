import param
import panel as pn
from matplotlib.colors import Colormap, Normalize
import proplot as pplt

class Opts(param.Parameterized):

    def __init__(self, **params):
        super().__init__(**params)
        self._excluded_from_opts = ['name']

        self.widgets = self.generate_widgets()

    @property
    def panel(self):
        return pn.Column(*self.widgets.values())

    @property
    def opts(self):
        opts = {name: self.param[name] for name in self.param if name not in self._excluded_from_opts}
        return opts

    def generate_widgets(self, **kwargs):
        """returns a dict with keys parameter names and values default mapped widgets"""
        #todo base class?

        names = [p for p in self.param if self.param[p].precedence is None or self.param[p].precedence > 1]
        widgets = pn.Param(self.param, show_name=False, show_labels=True, widgets=kwargs)

        return {k: v for k, v in zip(names[1:], widgets)}


class CmapOpts(Opts):

    cmap = param.ClassSelector(default=pplt.Colormap('viridis'), class_=Colormap)
    norm = param.ClassSelector(default=pplt.Norm('linear', 0., 1.), class_=Normalize)
    clim = param.Tuple((0., 1.), length=2)

    def __init__(self, **params):
        super().__init__(**params)
        self._excluded_from_opts += ['norm']
        self._norm_updated()

    @param.depends('norm', watch=True)
    def _norm_updated(self):
        self.clim = self.norm.vmin, self.norm.vmax
