import itertools
import warnings

import numpy as np
import pandas as pd
import panel as pn
import param
from param.parameterized import default_label_formatter

from pyhdx.support import autowrap, make_tuple
from pyhdx.web.sources import Source
from pyhdx.web.cache import Cache

# ABC
class Transform(param.Parameterized):
    """Gets data and applies transform"""

    _type = "base"

    widgets = param.Dict(default={})

    updated = param.Event()

    redrawn = param.Event(
        doc="event gets triggered when widgets are changed and the controller needs to redraw them"
    )

    _hash = param.Integer(doc="Hash of current transform state")

    _cache = param.ClassSelector(default=Cache(), class_=Cache)

    def __init__(self, **params):
        super().__init__(**params)

    # perhaps htey should all be private to prevent namespace collision with filter options
    @property
    def source_hash(self):
        return self.source.hash

    @property
    def hash_key(self):
        """hashable key describing the transform"""
        excluded = ["updated", "source", "widgets"]
        return tuple(
            (item, make_tuple(val))
            for item, val in self.param.get_param_values()
            if not (item.startswith("_") or item in excluded)
        )

    @property
    def hash(self):
        tup = (*self.hash_key, self.source_hash)

        return hash(tup)

    def update_hash(self):
        if self.hash == self._hash:
            return False
        else:
            self._hash = self.hash
            return True

    def update(self):
        if self.update_hash():
            self._update_options()  # todo this shouldnt be here
            self.updated = True


class MultiTransform(Transform):
    """transform which takes multiple sources as input"""

    _type = None

    sources = param.Dict(doc="Dict of sources the transform takes as input")

    @property
    def source_hash(self):
        return tuple(source.hash for source in self.sources.values())


class SelectTransform(MultiTransform):

    _type = "select"

    value = param.Selector()

    labels = param.List(
        default=[],
        doc="Optional list of labels for sources to display as options in selectors",
    )

    def __init__(self, **params):
        super().__init__(**params)
        for k, src in self.sources.items():
            src.param.watch(self.update, ["updated"])

        self.labels = self.labels or list(self.sources.keys())
        self.param["value"].objects = self.labels
        if not self.value:
            self.value = self.labels[0]
        self.redraw()

    @property
    def source_hash(self):
        return self.sources[self.value].hash

    # todo automagic widgets from labels
    def _make_widgets(self):
        return {"value": pn.Param(self.param, parameters=["value"], show_name=False)[0]}

    def redraw(self):
        self.widgets = self._make_widgets()

    def get(self):
        return self.sources[self.value].get()

    @pn.depends("value", watch=True)
    def update(self, *events):
        if self.update_hash():
            self.updated = True


class TableSourceTransform(Transform):
    """transform which picks the correct table from the source"""

    _type = "table_source"

    source = param.ClassSelector(class_=Source)

    table = param.Selector(
        default=None,
        doc="""
      The table being transformed. """,
    )

    def __init__(self, table_options=None, **params):
        self.table_options = table_options
        super().__init__(**params)
        self.widgets = {"table": pn.pane.panel(self.param.table)}

        if self.table_options:
            self._update_options()

    # todo allow auto generate widgets as in control panels /  views

    def get(self):
        df = self.source.get_table(
            self.table
        )  # returns None on KeyError #todo change to source.get_table
        return df

    @property
    def source_hash(self):
        return self.source.hashes.get(self.table, hash(None))

    def _update_options(self):
        options = self.source.get_tables()
        if self.table_options:
            options = [t for t in options if t in self.table_options]
        self.param["table"].objects = options
        if not self.table and options:
            self.table = options[0]

    @param.depends("source.updated", "table", watch=True)
    def update(self):
        self._update_options()
        if self.update_hash():
            self.updated = True


class AppTransform(Transform):
    """transform which acts on previous transforms in a chain. source is also a transform"""

    _type = None  # None type cannot be used in apps directly

    source = param.ClassSelector(class_=Transform)

    def transform(self):
        """get source data, apply transform, return result"""
        return self.source.get()

    def get(self):
        """method called to get the dataframe"""
        if self.hash in self._cache:
            return self._cache[self.hash]
        else:
            data = self.transform()
            self._cache[self.hash] = data
            return data

    @param.depends("source.updated", watch=True)
    def update(self):
        self.updated = True


class CrossSectionTransform(AppTransform):

    _type = "cross_section"

    pd_function = "xs"

    key = param.List()

    axis = param.Integer(1, bounds=[0, 1])

    n_levels = param.Integer(
        0, doc="Number of levels. negative to count from the back of the list"
    )

    level = param.List()

    drop_level = param.Boolean(
        True
    )  # note that there is a probable? bug when the xs result has one column not all levels are dropped

    names = param.List(None, doc="List of label names for widgets")

    empty_select = param.Boolean(
        default=False,
        doc="""
        Add an option to Select widgets to indicate select all on this level.""",
    )

    def __init__(self, **params):
        super().__init__(**params)
        self.index = (
            None  # index is the df index which determines the selector's options
        )
        self.update()

    @param.depends("source.updated", watch=True)
    def update(self):

        if self.update_hash():
            # todo remove watchers when new transforms are created?
            old_index = self.index
            df = self.source.get()

            if df is None:
                return
            self.index = df.columns if self.axis else df.index
            self._names = self.names or self.index.names

            if old_index is not None and self.index.nlevels == old_index.nlevels:
                # no redraw needed, only update selectors options
                options = list(self.index.unique(level=0))
                self.selectors[0].options = options
                self.selectors[0].param.trigger("value")
                for name, selector in zip(self._names, self.selectors):
                    selector.name = name  # todo requires testing if the names are really updated or not (they arent)
                    selector.label = name  # todo requires testing if the names are really updated or not
                    self.redrawn = True
            else:
                self.redraw()

            self.updated = True

    def redraw(self):
        # create new widgets
        if self.n_levels <= 0:
            n_levels = len(self._names) + self.n_levels
        else:
            n_levels = self.n_levels

        self.widgets = {
            name: pn.widgets.Select(name=default_label_formatter(name))
            for name in self._names[:n_levels]
        }

        # todo perhaps do self.param.add_parameter?
        self.selectors = list(self.widgets.values())
        for selector in self.selectors:
            selector.param.watch(self._selector_changed, ["value"], onlychanged=True)

        options = list(self.index.get_level_values(0).unique())
        if self.empty_select:  # todo use Nonetype? -> allow_none kwarg for Select?
            options = ["None"] + options
        self.selectors[0].options = options
        self.selectors[0].value = options[0]

        self.redrawn = True

    def transform(self):
        df = self.source.get()
        if df is None:
            return df

        kwargs = self.pd_kwargs
        # drop level bugged? https://github.com/pandas-dev/pandas/issues/6507
        df = df.xs(**kwargs)
        return df

    def _selector_changed(self, *events):
        # this sends multiple updated events as it triggers changes in other selectors
        for event in events:
            try:
                current_index = self.selectors.index(
                    event.obj
                )  # Index of the selector which was changed
            except ValueError:
                # I suspect this happens when there is redraw and the selector which triggered the change is no longer
                # in existence. So we can probaly ignore the event
                return

            try:  # try/except when we are at the last selector
                next_selector = self.selectors[current_index + 1]

                # Determine key/level to obtain new index which gives options for the next slider
                values = [
                    selector.value for selector in self.selectors[: current_index + 1]
                ]
                key = [value if value != "None" else slice(None) for value in values]
                level = list(range(current_index + 1))
                bools, current_columns = self.index.get_loc_level(key=key, level=level)
                options = list(
                    current_columns.get_level_values(0).unique()
                )  # current_columns.unique(level=0)
                if self.empty_select:
                    options = ["None"] + options
                next_selector.options = options
                if (
                    next_selector.value is None
                ):  # If the selector was not set yet, set it to 'None'
                    next_selector.value = options[0]
            except IndexError:
                pass

        # set the df
        all_values = [selector.value for selector in self.selectors]
        self.key = [value if value != "None" else slice(None) for value in all_values]
        self.level = list(range(len(all_values)))

        # TODO:
        # if self.update_hash() ?

        self.updated = True

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function"""
        return dict(
            key=tuple(self.key),
            axis=self.axis,
            level=self.level,
            drop_level=self.drop_level,
        )


class ApplyCmapOptTransform(AppTransform):

    _type = "apply_cmap_opt"

    opts = param.Selector(
        doc="cmap opts list to choose from", label="Color transform", objects=[],
    )  # todo refactor to cmap_opt? color transform?

    # @clasmethod
    # def check_args(... )  #todo method for constructor to see if the supplied kwargs are correct for this object

    def __init__(self, opts, **params):  # opts: list of opts objects
        warnings.warn(
            "ApplyCmapOptTransform does not implement hashing", NotImplementedError
        )
        self._opts_dict = {o.name: o for o in opts}
        opts = list(self._opts_dict.keys())
        params["opts"] = opts
        super().__init__(**params)
        self.widgets = {"opts": pn.pane.panel(self.param.opts)}

    def get(self):
        df = self.source.get()  # todo refactor df to data as it can also be a series?
        if df is None:
            return None

        if df.columns.size != 1:
            raise ValueError("Invalid number of columns, must be 1")

        pd_series = df[df.columns[0]]
        opts_obj = self._opts_dict[self.opts]
        if pd.api.types.is_numeric_dtype(pd_series):
            colors = opts_obj.apply(pd_series)

            return colors
        else:
            return None

    @param.depends("opts", watch=True)
    def _opts_changed(self):
        self.updated = True  # opts options are the same but selection changed, signal

    @param.depends("source.updated", watch=True)
    def update(self):
        pd_series = self.source.get()
        # todo just show all, later deal with setting the correct one? (infer from previous transform setting)
        if pd_series is None:
            # with param.parameterized.discard_events(self):
            self.param["opts"].objects = []
            self.opts = None
        else:
            options = list(
                self._opts_dict.keys()
            )  # this needs updating as opts_dict is static
            # with turn off param triggers, then update (unexpected results)
            # with param.parameterized.discard_events(self):
            self.param[
                "opts"
            ].objects = options  # TODO: pr/issue:? when setting objects which does not include the current setting selector is not reset?
            if self.opts is None and options:
                self.opts = options[0]
            elif self.opts not in options and options:  # todo or
                self.opts = options[0]
            elif not options:
                self.opts = None

        self.updated = True


class RectangleLayoutTransform(AppTransform):
    """
    Takes data with peptide start, end positions and transforms it to bottom, left, right, top positions of rectangles.

    By default, for a given interval start, end (inclusive, exclusive) a rectangle is returned with a height of spanning from
    (start - 0.5) to (end - 0.5)

    """

    _type = "rectangle_layout"

    left = param.String("start", doc="Field name to use for left coordinate")
    right = param.String("end", doc="Field name to use for for the right coordinate")

    height = param.Integer(1, doc="Height of the rectangles", constant=True)
    value = param.String(
        None, doc="Optional field name to pass through as value column"
    )
    passthrough = param.List(
        [], doc="Optional field names to pass through (for hovertools)"
    )

    wrap = param.Integer(
        doc="Amount of peptides to plot on the y axis before wrapping around to y axis top"
    )
    step = param.Integer(
        5,
        bounds=(1, None),
        doc="Step size used for finding 'wrap' when its not specified",
    )
    margin = param.Integer(
        4, doc="Margin space to keep between peptides when finding 'wrap'"
    )

    def transform(self):
        df = self.source.get()
        if df is None:
            return None

        df = df.dropna(subset=[self.left, self.right]).sort_values(
            by=[self.left, self.right]
        )

        # todo fix types to be ints in the df
        left = df[self.left].to_numpy(dtype=int)
        right = df[self.right].to_numpy(dtype=int)

        if not self.wrap:
            wrap = autowrap(left, right, margin=self.margin, step=self.step)
        else:
            wrap = self.wrap

        # Order of columns determines their role, not the names
        columns = [
            "x0",
            "y0",
            "x1",
            "y1",
        ]  # bottom-left (x0, y0) and top right (x1, y1)
        output_df = pd.DataFrame(index=df.index, columns=columns)
        output_df["x0"] = left - 0.5
        output_df["x1"] = right - 0.5
        cycle = itertools.cycle(
            range(self.height * wrap, 0, -self.height)
        )  # Starts at y top, cyles with wrap
        yvals = np.array(
            list(itertools.islice(cycle, len(df)))
        )  # Repeat cycle until the correct length
        output_df["y0"] = yvals - self.height
        output_df["y1"] = yvals

        if self.passthrough:
            for item in self.passthrough:
                assert item not in [
                    "value",
                    "index",
                ], "Invalid field name, 'index' and 'value' names are reserved"
                output_df[item] = df[item]
        output_df["index"] = df.index

        return output_df


class GenericTransform(AppTransform):
    _type = "generic"

    pd_function = param.String()

    kwargs = param.Dict(doc="dict of additional kwargs")

    def __init__(self, **params):
        self.kwargs = {k: v for k, v in params.items() if k not in self.param}
        super().__init__(**{k: v for k, v in params.items() if k in self.param})

    def transform(self):
        df = self.source.get()
        if df is None:
            return df
        func = getattr(df, self.pd_function)
        df = func(**self.pd_kwargs)

        return df

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function"""
        return self.kwargs


class DropLevelTransform(GenericTransform):
    _type = "droplevel"

    pd_function = "droplevel"


class RenameTransform(GenericTransform):
    _type = "rename"

    pd_function = "rename"

    def __init__(self, **params):
        self.kwargs = {k: v for k, v in params.items() if k not in self.param}
        super().__init__(**params)

    def transform(self):
        df = self.source.get()
        if df is None:
            return None
        kwargs = self.pd_kwargs.copy()
        columns = kwargs.pop("columns", None)
        if isinstance(columns, list):
            columns = {old: new for old, new in zip(df.columns, columns)}
        # todo same for index/rows?
        func = getattr(df, self.pd_function)
        df = func(columns=columns, **kwargs)

        return df


class ResetIndexTransform(GenericTransform):
    _type = "reset_index"

    pd_function = "reset_index"


class RescaleTransform(GenericTransform):
    """Rescale a single column"""

    _type = "rescale"

    pd_function = "assign"

    columns = param.List(doc="Name of the columns to rescale")

    scale_factor = param.Number(1.0)

    def transform(self):  # todo perhaps some kind of decorator that returns nonealwasy?
        df = self.source.get()
        if df is None:
            return None
        for column in self.columns:
            df = df.assign(**{column: lambda x: x[column] * self.scale_factor})

        return df


class PivotTransform(GenericTransform):
    _type = "pivot"

    pd_function = "pivot"

    index = param.ClassSelector(class_=[str, list])

    columns = param.ClassSelector(class_=[str, list])

    values = param.ClassSelector(class_=[str, list])

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function"""
        # todo get_params func which finds the correct params here
        return dict(index=self.index, columns=self.columns, values=self.values)


class StackTransform(GenericTransform):
    _type = "stack"

    pd_function = "stack"

    level = param.Integer(
        -1
    )  # actually int, str, or list of these, default -1 (last level)
    # level = param.ClassSelector(_class=[int, str, list])  # where list is list of (str | int)

    dropna = param.Boolean(True)

    @property
    def pd_kwargs(self):
        """kwargs to pass to pandas function"""
        # todo get_params func which finds the correct params here
        return dict(level=self.level, dropna=self.dropna)


class ConcatTransform(AppTransform):
    """app transform which combines multiple dfs"""

    def __init__(self, **params):
        raise NotImplementedError()


# resampletransform?
class SampleTransform(AppTransform):
    """subsamples dataframe along specified axis"""

    _type = "sample"

    n = param.Integer(
        None,
        bounds=(0, None),
        doc="Number of subsamples to return. Incompatible with frac",
    )

    frac = param.Number(
        None,
        bounds=(0.0, 1.0),
        doc="Fraction of subsamples to return. Incompatible with n",
    )

    random = param.Boolean(
        False, doc="whether to sample randomly or equidistanly spaced"
    )

    axis = param.Number(0, inclusive_bounds=(0, 1))

    def transform(self):
        df = self.source.get()
        if df is None:
            return None

        if (self.frac is None) != (self.n is not None):
            raise ValueError("Must specify either 'n' or 'frac'")

        if self.random:
            df = df.sample(
                n=self.n, frac=self.frac, axis=self.axis
            )  # todo allow other kwargs?
        else:
            size = df.shape[self.axis]
            if self.n is not None:
                step = size // self.n
            else:
                step = int(1 / self.frac)

            step = max(step, 1)
            slice_step = slice(None, None, step)
            slice_all = slice(None)

            col_slicer = slice_step if self.axis else slice_all
            row_slicer = slice_all if self.axis else slice_step

            df = df.iloc[row_slicer, col_slicer]

        return df


class PipeTransform(AppTransform):
    """applies a list of pandas functions

    (not exactly the same as df.pipe)

    """

    _type = "pipe"

    pipe = param.List()  # list of dicts

    def transform(self):
        df = self.source.get()
        if df is None:
            return None
        for d in self.pipe:
            func = getattr(df, d["function"])
            args = d.get("args", [])
            kwargs = d.get("kwargs", {})

            df = func(*args, **kwargs)

        return df
