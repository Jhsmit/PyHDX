from copy import copy

import param
import panel as pn
from matplotlib.colors import Colormap, Normalize
import proplot as pplt

from pyhdx.plot import default_cmap_norm
from pyhdx.support import apply_cmap


class Opts(param.Parameterized):

    updated = param.Event()

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


class GenericOpts(Opts):

    def __init__(self, **params):
        self.kwargs = {k: v for k, v in params.items() if k not in self.param}
        super().__init__(**{k: v for k, v in params.items() if k in self.param})

    @property
    def opts(self):
        return self.kwargs



class CmapOpts(Opts):

    cmap = param.ClassSelector(default=None, class_=Colormap)

    norm = param.ClassSelector(default=None, class_=Normalize)
    # the stored norm here is the scaled one
    # scale factor is applied to apply norm to rescaled data

    clim = param.Tuple((0., 1.), length=2)

    sclf = param.Number(1., doc='scaling factor to apply')

    field = param.String(doc="optional field on which cmap works")

    def __init__(self, **params):
        cmap = params.pop('cmap', None)
        cmap = pplt.Colormap(cmap) if cmap else cmap
        params['cmap'] = cmap
        super().__init__(**params)
        self._excluded_from_opts += ['norm', 'sclf']  # perhaps use leading underscore to exclude?

        if self.cmap is None and self.norm is None and self.field is not None:
            self.cmap, self.norm = default_cmap_norm(self.field)
        elif self.field is None:
            self.cmap = pplt.Colormap('viridis')
            self.norm = pplt.Norm('linear', 0., 1.)
        self._norm_updated()

    @property
    def opts(self):
        names = ['cmap', 'clim']
        opts = {name: self.param[name] for name in names}
        return opts

    @param.depends('norm', watch=True)
    def _norm_updated(self):
        self.clim = self.norm.vmin, self.norm.vmax

    def apply(self, data):
        """apply cmap / norm to data (pd series or df)"""
        norm = copy(self.norm)
        norm.vmin *= self.sclf
        norm.vmax *= self.sclf
        return apply_cmap(data, self.cmap, norm)

    @param.depends('norm', 'cmap', watch=True)
    def update(self):
        self.updated = True