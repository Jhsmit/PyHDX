from copy import copy
from functools import partial, reduce

import panel as pn
import param
import proplot as pplt
from matplotlib.colors import Colormap, Normalize

from pyhdx.plot import CMAP_NORM_DEFAULTS
from pyhdx.support import apply_cmap


# todo baseclass widget generating thingy
class OptsBase(param.Parameterized):

    _type = None

    updated = param.Event()

    def __init__(self, **params):
        super().__init__(**params)
        self._excluded_from_opts = ["name"]  # todo remove this and opts property

        self.widgets = self.generate_widgets()

    @property
    def panel(self):
        return pn.Column(*self.widgets.values())

    @property
    def opts(self):
        opts = {
            name: self.param[name]
            for name in self.param
            if name not in self._excluded_from_opts
        }
        return opts

    def generate_widgets(self, **kwargs):
        """returns a dict with keys parameter names and values default mapped widgets"""
        # todo base class?

        names = [
            p
            for p in self.param
            if self.param[p].precedence is None or self.param[p].precedence > 1
        ]
        widgets = pn.Param(
            self.param, show_name=False, show_labels=True, widgets=kwargs
        )

        return {k: v for k, v in zip(names[1:], widgets)}


class GenericOpts(OptsBase):

    _type = "generic"

    hooks = param.List()

    def __init__(self, **params):
        self.kwargs = {k: v for k, v in params.items() if k not in self.param}
        super().__init__(**{k: v for k, v in params.items() if k in self.param})

    def hooks_factory(self):
        def hook(hooks, plot, element):
            for hook_spec in hooks:
                handle = plot.handles[hook_spec["handle"]]
                rsetattr(handle, hook_spec["attr"], hook_spec["value"])

        f = partial(hook, self.hooks)

        return f

    @property
    def opts(self):
        #  self.kwargs.update({'hooks': [self.hooks_factory()]})
        return {"hooks": [self.hooks_factory()], **self._parse_kwargs(self.kwargs)}

    @staticmethod
    def _parse_kwargs(kwargs):
        out = {}
        for k, v in kwargs.items():
            if k in ["xlim", "ylim"] and isinstance(v, list):
                out[k] = tuple(v)
            elif k in ["padding"] and isinstance(v, list):
                out[k] = to_tuple(v)
            else:
                out[k] = v

        return out


# https://stackoverflow.com/questions/27049998/convert-a-mixed-nested-list-to-a-nested-tuple/27050037#27050037
def to_tuple(lst):
    return tuple(to_tuple(i) if isinstance(i, list) else i for i in lst)


# https://stackoverflow.com/questions/31174295/getattr-and-setattr-on-nested-subobjects-chained-properties
def rsetattr(obj, attr, val):
    pre, _, post = attr.rpartition(".")
    return setattr(rgetattr(obj, pre) if pre else obj, post, val)


def rgetattr(obj, attr, *args):
    def _getattr(obj, attr):
        return getattr(obj, attr, *args)

    return reduce(_getattr, [obj] + attr.split("."))


class HooksOpts(OptsBase):

    _type = "hooks"

    hooks = param.List()

    def hooks_factory(self):
        def hook(hooks, plot, element):
            for hook_spec in hooks:
                handle = plot.handles[hook_spec["handle"]]
                rsetattr(handle, hook_spec["attr"], hook_spec["value"])

        f = partial(hook, self.hooks)

        return f

    @property
    def opts(self):
        opts = {"hooks": [self.hooks_factory()]}
        return opts


class CmapOpts(OptsBase):

    _type = "cmap"

    cmap = param.ClassSelector(default=None, class_=Colormap)

    norm = param.ClassSelector(default=None, class_=Normalize)
    # the stored norm here is the unscaled one
    # scale factor is applied to clim and norm_scaled

    clim = param.Tuple((0.0, 1.0), length=2)

    sclf = param.Number(1.0, doc="scaling factor to apply")  # curent: 1e+3

    field = param.String(doc="field on which cmap works")

    def __init__(self, rename=True, invert=False, **params):
        # todo from_spec constructor method for this kind of logic
        cmap = params.pop("cmap", None)
        cmap = pplt.Colormap(cmap) if cmap else cmap
        params["cmap"] = cmap
        super().__init__(**params)
        self._excluded_from_opts += [
            "norm",
            "sclf",
        ]  # perhaps use leading underscore to exclude?

        if self.cmap is None and self.norm is None and self.field is not None:
            cmap, norm = CMAP_NORM_DEFAULTS[self.field]
        elif self.field is None:
            cmap = pplt.Colormap("viridis")
            norm = pplt.Norm("linear", 0.0, 1.0)

        self.norm = norm
        self._cmap = cmap  # unreversed cmap

        if rename:
            cmap.name = self.field + "_default"
        if invert:
            cmap = cmap.reversed()

        self.cmap = cmap

        # self._norm_updated()

        # self.cmap = self.cmap.reversed()

    @property
    def opts(self):
        names = ["cmap", "clim"]
        opts = {name: self.param[name] for name in names}
        return opts

    @property
    def norm_scaled(self):
        norm = copy(self.norm)
        norm.vmin *= self.sclf
        norm.vmax *= self.sclf

        return norm

    @norm_scaled.setter
    def norm_scaled(self, norm):
        _norm = copy(norm)
        _norm.vmin /= self.sclf
        _norm.vmax /= self.sclf

        self.norm = _norm

    @param.depends("norm", watch=True)
    def _norm_updated(self):
        self.clim = self.norm.vmin * self.sclf, self.norm.vmax * self.sclf
        # todo invert bool?
        # self.clim = self.norm.vmax*self.sclf, self.norm.vmin*self.sclf,

    def apply(self, data):
        """apply cmap / norm to data (pd series or df)"""
        # norm = copy(self.norm)
        # norm.vmin *= self.sclf
        # norm.vmax *= self.sclf
        return apply_cmap(data, self.cmap, self.norm)

    @param.depends("norm", "cmap", watch=True)
    def update(self):
        self.updated = True
