import itertools
import logging
import time
from functools import partial
from itertools import groupby, count

import holoviews as hv
import numpy as np
import pandas as pd
import panel as pn
import param
from holoviews.streams import Pipe
from hvplot import hvPlotTabular
from panel.pane.base import PaneBase

from pyhdx.support import hex_to_rgb
from pyhdx.web.pane import PDBeMolStar, REPRESENTATIONS
from pyhdx.web.sources import Source
from pyhdx.web.transforms import Transform
from pyhdx.web.widgets import LoggingMarkdown, COLOR_SCHEMES, NGL
from pyhdx.web.widgets import REPRESENTATIONS as NGL_REPRESENTATIONS
from pyhdx.web.opts import CmapOpts


class View(param.Parameterized):
    """Base view object.

    Inspired by Holoviz Lumen's View objects"""

    _type = None

    opts = param.List(
        default=[], doc="list of opts dicts to apply on the plot", precedence=-1
    )

    dependencies = param.List(
        default=[],
        precedence=-1,
        doc="Additional dependencies which trigger update when their `updated` event fires",
    )

    def __init__(self, **params):
        super().__init__(**params)
        # todo allow for kwargs to be passed to DynamicMap's func

        for dep in self.dependencies:
            dep.param.watch(self._update, ["updated"])

        self.widgets = self.make_dict()

        self._panel = None
        self._updates = None  # what does this do?

    def _update(self, *events):
        self.update()  # todo or just catch events in the update function?  (probably this)

    def make_dict(self):
        # todo baseclass filter/controller/view with the same widget generation logic?
        """dict of widgets to be shown
        override this method to get custom mapping

        """
        return self.generate_widgets()

    def generate_widgets(self, **kwargs):
        """returns a dict with keys parameter names and values default mapped widgets"""

        names = [
            p
            for p in self.param
            if self.param[p].precedence is None or self.param[p].precedence > 1
        ]
        widgets = pn.Param(
            self.param, show_name=False, show_labels=True, widgets=kwargs
        )

        return {k: v for k, v in zip(names[1:], widgets)}

    def get_data(self):  # refactor get?
        """
        Queries the Source

        Returns
        -------
        DataFrame
            The queried table after filtering and transformations are
            applied.
        """
        # if self._cache is not None:
        #     return self._cache

        df = self.source.get()

        return df

    def _update_panel(self, *events):
        """
        Updates the cached Panel object and returns a boolean value
        indicating whether a rerender is required.
        """

        # todo clenaup
        if self._panel is not None:
            self._cleanup()
            self._updates = self._get_params()
            if self._updates is not None:
                return False
        self._panel = self.get_panel()
        return True

    @property  # todo move to init?
    def opts_dict(self):
        # Combine all opts and merge overlapping lists of opts
        # (currently to merge hooks, might be unwanted behaviour for others)
        opts_dict = {}
        for d in self.opts:  # Iterate over Opts on this view
            for k, v in d.opts.items():  # Iterate over all dict items in the Opt
                if k in opts_dict:
                    if isinstance(v, list) and isinstance(opts_dict[k], list):
                        combined_list = opts_dict[k] + v
                        opts_dict[k] = combined_list
                    else:
                        raise ValueError(
                            f"Overlapping key {k!r} in opt {d.name!r} on view {self.name!r}"
                        )
                else:
                    opts_dict[k] = v

        return opts_dict


class hvView(View):

    source = param.ClassSelector(
        class_=(Source, Transform),
        constant=True,
        precedence=-1,
        doc="""
        The Source to query for the data.""",
    )

    _type = None

    def __init__(self, **params):
        super().__init__(**params)
        self._stream = None

    @param.depends("source.updated", watch=True)
    def update(self):
        """
        Triggers an update in the View.

        Returns
        -------
        stale : bool
            Whether the panel on the View is stale and needs to be
            rerendered.
        """
        # Skip events triggered by a parameter change on this View
        # own_parameters = [self.param[p] for p in self.param]
        # own_events = events and all(
        #     isinstance(e.obj, ParamFilter) and
        #     (e.obj.parameter in own_parameters or
        #     e.new is self._ls.selection_expr)
        #     for e in events
        # )
        # if own_events:
        #     return False
        # if invalidate_cache:
        #     self._cache = None
        if self._stream is None:
            return self._update_panel()
        if self.get_data() is not None:
            self._stream.send(self.get_data())
        return False

    def get_panel(self):
        kwargs = self._get_params()
        # interactive? https://github.com/holoviz/panel/issues/1824
        return pn.pane.HoloViews(linked_axes=False, **kwargs)  # linked_axes=False??

    def _get_params(self):
        df = self.get_data()
        if df is None:
            df = self.empty_df

        self._stream = Pipe(data=df)
        return dict(
            object=self.get_plot(), sizing_mode="stretch_both"
        )  # todo update sizing mode

    @property
    def panel(self):
        if isinstance(self._panel, PaneBase):
            pane = self._panel
            if len(pane.layout) == 1 and pane._unpack:
                return pane.layout[0]
            return pane._layout
        return self._panel


class hvPlotView(hvView):
    _type = "hvplot"

    kind = param.String()

    def __init__(self, **params):
        self.kwargs = {k: v for k, v in params.items() if k not in self.param}
        super().__init__(**{k: v for k, v in params.items() if k in self.param})
        self._stream = None

    def get_plot(self):
        """

        Parameters
        ----------
        df

        Returns
        -------

        """

        def func(data, kind, **kwargs):
            return hvPlotTabular(data)(kind=kind, **kwargs)

        pfunc = partial(func, kind=self.kind, **self.kwargs)

        plot = hv.DynamicMap(pfunc, streams=[self._stream])
        plot = plot.apply.opts(**self.opts_dict)

        return plot

    @property
    def empty_df(self):
        df = pd.DataFrame({"null": [np.nan], "y2": [np.nan]})
        return df


class hvCurveView(hvView):
    _type = "curve"

    x = param.String(
        None, doc="The column to render on the x-axis."
    )  # todo these should be selectors

    y = param.String(None, doc="The column to render on the y-axis.")

    def __init__(self, **params):
        # baseclass??
        self._stream = None
        super().__init__(**params)

    def get_plot(self):
        """

        Parameters
        ----------
        df

        Returns
        -------

        """

        func = partial(hv.Curve, kdims=self.kdims, vdims=self.vdims)
        plot = hv.DynamicMap(func, streams=[self._stream])
        plot = plot.apply.opts(**self.opts_dict)

        return plot

    @property
    def kdims(self):
        return [self.x] if self.x is not None else None

    @property
    def vdims(self):
        return [self.y] if self.y is not None else None

    @property
    def empty_df(self):
        columns = (self.kdims or ["x"]) + (self.vdims or ["y"])
        return pd.DataFrame([[np.nan] * len(columns)], columns=columns)


class hvScatterAppView(hvView):
    _type = "scatter"

    x = param.String(
        None, doc="The column to render on the x-axis."
    )  # todo these should be selectors

    y = param.String(None, doc="The column to render on the y-axis.")

    def __init__(self, **params):
        self._stream = None
        super().__init__(**params)

    def get_plot(self):
        """

        Parameters
        ----------
        df

        Returns
        -------

        """

        func = partial(hv.Scatter, kdims=self.kdims, vdims=self.vdims)
        plot = hv.DynamicMap(func, streams=[self._stream])
        plot = plot.apply.opts(**self.opts_dict)

        return plot

    @property
    def kdims(self):
        return [self.x] if self.x is not None else None

    @property
    def vdims(self):
        return [self.y] if self.y is not None else None

    @property
    def empty_df(self):
        dic = {self.x or "x": [], self.y or "y": []}
        if "color" in self.opts_dict:
            dic[self.opts_dict["color"]] = []
        return pd.DataFrame(dic)


class hvRectanglesAppView(hvView):
    _type = "rectangles"

    x0 = param.String("x0")

    x1 = param.String("x1")

    y0 = param.String("y0")

    y1 = param.String("y1")

    vdims = param.List(["value"])

    def __init__(self, **params):
        # todo left and right cannot be none?
        super().__init__(**params)

    def get_plot(self):
        """
        Dataframe df must have columns x0, y0, x1, y1 (in this order) for coordinates
        bottom-left (x0, y0) and top right (x1, y1). Optionally a fifth value-column can be provided for colors

        Parameters
        ----------
        df

        Returns
        -------

        """

        func = partial(hv.Rectangles, kdims=self.kdims, vdims=self.vdims)
        plot = hv.DynamicMap(func, streams=[self._stream])

        if self.opts_dict:
            plot = plot.apply.opts(**self.opts_dict)

        return plot

    @property
    def kdims(self):
        return [self.x0, self.y0, self.x1, self.y1]

    @property
    def empty_df(self):
        columns = self.kdims + self.vdims
        return pd.DataFrame([[0] * len(columns)], columns=columns)


class hvErrorBarsAppView(hvView):
    _type = "errorbars"

    pos = param.String(
        "x", doc="Positions of the errobars, x-values for vertical errorbars"
    )

    value = param.String(
        "y", doc="Values of the samples, y-values for vertical errorbars"
    )

    err = param.String(None, doc="Error values in both directions")

    err_pos = param.String(None, doc="Error values in positive direction")

    err_neg = param.String(None, doc="Error values in negative direction")

    horizontal = param.Boolean(False, doc="error bar direction")

    def __init__(self, **params):

        # todo left and right cannot be none?
        super().__init__(**params)

    def get_plot(self):
        """
        Dataframe df must have columns x0, y0, x1, y1 (in this order) for coordinates
        bottom-left (x0, y0) and top right (x1, y1). Optionally a fifth value-column can be provided for colors

        Parameters
        ----------
        df

        Returns
        -------

        """

        func = partial(
            hv.ErrorBars, kdims=self.kdims, vdims=self.vdims, horizontal=self.horizontal
        )
        plot = hv.DynamicMap(func, streams=[self._stream])

        if self.opts_dict:
            plot = plot.apply.opts(**self.opts_dict)

        return plot

    @property
    def vdims(self):
        if self.err is not None and self.err_pos is None and self.err_neg is None:
            return [self.value, self.err]
        elif self.err is None and self.err_pos is not None and self.err_neg is not None:
            return [self.value, self.err_pos, self.err_neg]
        else:
            raise ValueError(
                "Must set either only 'err' or both 'err_pos' and 'err_neg'"
            )

    @property
    def kdims(self):
        return [self.pos]

    @property
    def empty_df(self):
        columns = self.kdims + self.vdims
        return pd.DataFrame([[0] * len(columns)], columns=columns)


class hvOverlayView(View):
    _type = "overlay"

    views = param.List(doc="List of view instances to make overlay")

    def update(self):
        self._update_panel()

    def _cleanup(self):
        pass

    def get_plot(self):
        items = [view.get_plot() for view in self.views]
        plot = hv.Overlay(
            items
        ).collate()  # todo is collate always needed? Does it always return a DynamicMap? -> generalize

        if self.opts_dict:
            plot = plot.apply.opts(**self.opts_dict)

        return plot

    def _get_params(self):
        return dict(object=self.get_plot(), sizing_mode="stretch_both")

    def get_panel(self):
        kwargs = self._get_params()
        # interactive? https://github.com/holoviz/panel/issues/1824
        return pn.pane.HoloViews(**kwargs)

    @property
    def panel(self):
        if isinstance(self._panel, PaneBase):
            pane = self._panel
            if len(pane.layout) == 1 and pane._unpack:
                return pane.layout[0]
            return pane._layout
        return self._panel


class PDBeMolStarColorView(View):
    _type = "pdbemolstar_colors"

    sources = param.Dict(
        doc="Dict of sources for this view. "
        "should be: pdb: PDBSource, color: TableSource (single-column tables)"
    )

    visual_style = param.Selector(default="cartoon", objects=REPRESENTATIONS)

    lighting = param.Selector(
        default="matte",
        objects=["flat", "matte", "glossy", "metallic", "plastic"],
        doc="Set the lighting",
    )

    spin = param.Boolean(default=False, doc="Toggle object spin")

    highlight_color = param.Color(
        default="#ff6699", doc="Color for mouseover highlighting"
    )

    reset = param.Action(lambda self: self._action_reset())

    def __init__(self, **params):
        # also accepts any PDBeMolstar kwargs
        # todo should generate widgets which can be displayed in the controller
        # get kwargs
        self.pdbe_kwargs = {k: v for k, v in params.items() if k not in self.param}
        super().__init__(**{k: v for k, v in params.items() if k in self.param})

        self.pdbe = None
        self._panel = pn.Row(sizing_mode="stretch_both")

        self.sources["pdb"].param.watch(self._pdb_updated, "updated")
        self.sources["color"].param.watch(self._color_updated, "updated")

        # field: opts for all cmap otps
        self._cmap_opts = {
            opt.field: opt for opt in self.opts if isinstance(opt, CmapOpts)
        }

    def _action_reset(self):
        data = {"camera": True}
        self.pdbe.reset(data)
        self._color_updated("event")  # todo new function without event arg

    # todo perhaps wrappertje @not_none('pdbe') which checks for None and returns
    def _update_params(self, *events):
        if self.pdbe is None:
            return
        rerender = ["visual_style", "lighting"]
        for event in events:
            if event.name in rerender:
                # color_data_dict = self.get_color_data()
                # color_data_dict["nonSelectedColor"] = color_data_dict.pop("non_selected_color")
                # self.pdbe.color_on_load = color_data_dict

                time.sleep(2)  # in the future bind to load complete event
                self._color_updated("event")  # update?

            setattr(self.pdbe, event.name, event.new)

    def get_panel(self):
        return self._panel

    def _cleanup(self):
        return None

    def _get_params(self):
        return None

    def init_pdbe(self, pdb_url):
        self.pdbe = PDBeMolStar(
            sizing_mode="stretch_both",
            custom_data={"url": pdb_url, "format": "pdb"},
            **self.pdbe_kwargs,
        )

        # link view params to pdbe params
        params = self.param.params().keys() & self.pdbe.param.params().keys() - {"name"}
        self.param.watch(self._update_params, list(params))
        self._panel[:] = [self.pdbe]

    def _pdb_updated(self, *events):
        pdb_url = self.sources["pdb"].get()

        if self.pdbe is None:
            self.init_pdbe(pdb_url)
        else:
            self.pdbe.custom_data = {"url": pdb_url, "format": "pdb"}

        # color_data_dict = self.get_color_data()
        # color_data_dict["nonSelectedColor"] = color_data_dict.pop("non_selected_color")
        # self.pdbe.color_on_load = color_data_dict

        time.sleep(2)  # todo link to done event
        self._color_updated("event")

    def get_color_data(self):
        df = self.sources["color"].get()
        if df is None:
            return None
        # there should be one column which matches one of the keys in the cmap otps dict
        matching_columns = set(df.columns) & self._cmap_opts.keys()
        if not matching_columns:
            # todo logging.getlogger etc etc
            print("No matching color opts were found")
            return None

        qty = matching_columns.pop()
        opts = self._cmap_opts[qty]
        r, g, b, a = opts.cmap.get_bad() * 255
        no_coverage = {"r": r, "g": g, "b": b}
        # pd.Series with colors, take entries with residue number index >= 1
        color_series = opts.apply(df[qty]).loc[1:]

        # Group subsequent residues with identical color values
        # Use these grouped value to generate data dictionary to pass to PDBeMolstar
        # to apply color scheme
        colors = color_series.values
        r_numbers = color_series.index.values
        data_list = []
        i = 0
        for key, grp in itertools.groupby(colors):
            size = sum(1 for x in grp)
            data_elem = {
                "start_residue_number": r_numbers[i],
                "end_residue_number": r_numbers[i + size - 1],
                "color": {k: v for k, v in zip("rgb", hex_to_rgb(key))},
            }
            data_list.append(data_elem)
            i += size

        return {"data": data_list, "non_selected_color": no_coverage}

    def _color_updated(self, *event):
        if self.pdbe is None:
            return

        color_data_dict = self.get_color_data()
        if color_data_dict is not None:
            self.pdbe.color(**color_data_dict)

        return self._update_panel()  # check how and why this is needed

    # this is called to initiate the view. perhaps should be removed / refacotred
    # its also triggered by any dependency trigger (in this case opts )
    # it is not triggered by sources triggering (which eg hvplot view does do)
    def update(self):
        self._color_updated()
        return self._update_panel()

    @property
    def panel(self):  # why the panebase unpack?
        if isinstance(self._panel, PaneBase):
            pane = self._panel
            if len(pane.layout) == 1 and pane._unpack:
                return pane.layout[0]
            return pane._layout
        return self._panel


class NGLColorView(View):
    _type = "ngl_colors"

    # todo additioal render options (fog etc)

    sources = param.Dict(
        doc="Dict of sources for this view. "
        "should be: pdb: PDBSource, color: TableSource (single-column tables)"
    )

    representation = param.Selector(default="cartoon", objects=NGL_REPRESENTATIONS)

    effect = param.Selector(
        default=None, objects=[None, "spin", "rock"], allow_None=True
    )

    color_scheme = param.Selector(default="custom", objects=COLOR_SCHEMES)

    custom_color_scheme = param.List(precedence=-1)

    background_color = param.Color(default="#F7F7F7")

    object = param.String("", doc="pdb string object", precedence=-1)

    def __init__(self, **params):
        # todo should generate widgets which can be displayed in the controller
        super(NGLColorView, self).__init__(**params)
        self._ngl = NGL(
            sizing_mode="stretch_both",  # todo sanitize order
            extension="pdb",
            effect=self.effect,
            color_scheme=self.color_scheme,
            representation=self.representation,
            object=self.object,
        )

        params = self.param.params().keys() & self._ngl.param.params().keys() - {"name"}
        self.param.watch(self._update_params, list(params))

        self.sources["pdb"].param.watch(self._pdb_updated, "updated")
        self.sources["color"].param.watch(self._color_updated, "updated")

        # field: opts for all cmap otps
        self._cmap_opts = {
            opt.field: opt for opt in self.opts if isinstance(opt, CmapOpts)
        }

    def _update_params(self, *events):
        for event in events:
            setattr(self._ngl, event.name, event.new)

    def get_panel(self):
        return self._ngl

    def _cleanup(self):
        return None

    def _get_params(self):
        return None

    def _pdb_updated(self, *events):
        pdb_string = self.sources["pdb"].get()
        self.object = pdb_string

    def _color_updated(self, *event):
        df = self.sources["color"].get()
        if df is None:
            return

        # there should be one column which matches one of the keys in the cmap otps dict
        matching_columns = set(df.columns) & self._cmap_opts.keys()
        if not matching_columns:
            # todo logging.getlogger etc etc
            print("No matching color opts were found")
            return

        qty = matching_columns.pop()
        opts = self._cmap_opts[qty]
        colors = opts.apply(df[qty])  # pd.Series with colors

        grp = colors.groupby(
            colors
        )  # Group to color and transform to ngl custom color scheme format
        color_list = []
        for c, pd_series in grp:
            result = [
                list(g)
                for _, g in groupby(
                    pd_series.index, key=lambda n, c=count(): n - next(c)
                )
            ]

            resi = " or ".join([f"{g[0]}-{g[-1]}" for g in result])
            color_list.append([c, resi])

        self.custom_color_scheme = color_list

        return self._update_panel()  # check how and why this is needed

    # this is called to initiate the view. perhaps should be removed / refacotred
    # its also triggerd by any dependency trigger (in this case opts)
    def update(self):
        self._color_updated()

        return self._update_panel()

    @property
    def panel(self):  # why the panebase unpack?
        if isinstance(self._panel, PaneBase):
            pane = self._panel
            if len(pane.layout) == 1 and pane._unpack:
                return pane.layout[0]
            return pane._layout
        return self._panel


class LoggingView(View):
    _type = "logging"

    logger = param.ClassSelector(
        logging.Logger, doc="Logger object to show in Log view"
    )

    level = param.Integer(
        default=10,
        doc="Logging level of the streamhandler redirecting logs to the view",
    )

    def __init__(self, *args, **params):
        super(LoggingView, self).__init__(**params)
        title = self.opts_dict.get("title", "Log Window")
        self.markdown = LoggingMarkdown(f"### {title} \n", sizing_mode="stretch_both")

        self.sh = logging.StreamHandler(self.markdown)
        self.sh.terminator = "  \n"
        self.sh.setLevel(self.level)
        formatter = logging.Formatter(
            "%(asctime)s [%(levelname)s]: %(message)s", "%Y-%m-%d %H:%M:%S"
        )
        self.sh.setFormatter(formatter)
        self.logger.addHandler(self.sh)

    @param.depends("level", watch=True)
    def _level_updated(self):
        self.sh.setLevel(self.level)

    @property
    def panel(self):
        return self.markdown

    def update(self, *events, invalidate_cache=True):
        pass
