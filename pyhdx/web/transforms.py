from lumen.transforms import Transform
import param
from matplotlib.colors import Colormap, Normalize
from pyhdx.support import rgb_to_hex, autowrap
import itertools
import pandas as pd
import numpy as np
import panel as pn


class WebAppTransform(param.Parameterized):  #todo subclass from Transform?

    updated = param.Event()

    def __init__(self, **params):
        super().__init__(**params)
        self.widgets = self.generate_widgets()

    @property
    def panel(self):
        return pn.Column(*self.widgets.values())
        # widget = self.widget.clone()
        # self.widget.link(widget, value='value', bidirectional=True)
        # return widget

    def generate_widgets(self, **kwargs):
        """returns a dict with keys parameter names and values default mapped widgets"""

        names = [p for p in self.param if self.param[p].precedence is None or self.param[p].precedence > 1]
        widgets = pn.Param(self.param, show_name=False, show_labels=True, widgets=kwargs)

        return {k: v for k, v in zip(names[1:], widgets)}


class RescaleTransform(Transform):
    """
    This transform takes a field from a table and transforms in with a scaling factors
    """

    transform_type = 'rescale'

    scale_factor = param.Number(doc='Scaling factor to multiply field values with')
    field = param.String(doc='Name of the field to transform')

    def apply(self, table):
        table = table.copy()  # Performance-wise this might not be ideal
        table[self.field] *= self.scale_factor

        return table


class AccumulateRegularizersTransform(Transform):
    """
    Very niche and temporary transform to accumulate reg losses to one column
    """

    transform_type = 'accumulate_regularizers'

    def apply(self, table):
        # first two columns are index and mse_loss?
        reg_total = table.iloc[:, 2:].sum(axis=1)
        reg_total.name = 'reg_loss'

        result = pd.concat([table.iloc[:, :2], reg_total], axis=1)

        return result


class ResetIndexTransform(Transform):

    level = param.ClassSelector(class_=(int, list, str), doc="""
        Only remove the given levels from the index. Removes all levels by default.""")

    drop = param.Boolean(default=False, doc="""
        Do not try to insert index into dataframe columns. This resets the index to the default integer index.""")

    col_level = param.ClassSelector(default=0, class_=(int, str), doc="""
        If the columns have multiple levels, determines which level the labels are inserted into. By default it is 
        inserted into the first level.""")

    col_fill = param.Parameter(default='', doc="""If the columns have multiple levels, determines how the other 
        levels are named. If None then the index name is repeated.""")

    transform_type = 'reset_index'

    def apply(self, table):
        return table.reset_index(level=self.level, drop=self.drop, col_level=self.col_level, col_fill=self.col_fill)


class SetIndexTransform(Transform):

    keys = param.Parameter(doc="""
        This parameter can be either a single column key, a single array of the same length as the calling DataFrame, 
        or a list containing an arbitrary combination of column keys and arrays. Here, “array” encompasses Series, 
        Index, np.ndarray, and instances of Iterator.""") ## label or array-like or list of labels/arrays

    drop = param.Boolean(default=True, doc="""
        Delete columns to be used as the new index.""")

    append = param.Boolean(default=False, doc="""
        Whether to append column to an existing index""")

    verify_integrity = param.Boolean(default=False, doc="""
        Check the new index for duplicates. Otherwise defer the check until necessary. Setting to False will improve 
        the performance of this method.""")

    transform_type = 'set_index'

    def apply(self, table):
        return table.set_index(self.keys, drop=self.drop, append=self.append, inplace=self.inplace,
                               verify_integrity=self.verify_integrity)


class ApplyCmapTransform(Transform, WebAppTransform):
    """
    This transform takes data from a specified field, applies a norm and color map, and adds the resulting colors
    in a new column
    """

    fields = param.Selector(doc='Fields to choose from to apply cmap to')

    field = param.String(doc='Name of the field to apply colors to')
    cmap = param.ClassSelector(Colormap)
    norm = param.ClassSelector(Normalize)
    color_column = param.String('color', doc='Name of the added color column')
    # default_color =

    transform_type = 'color'

    def __init__(self, **params):
        super().__init__(**params)

        #temporariy
        self.param['fields'].objects = ['deltaG', 'pfact']
        self.fields = 'deltaG'

    def apply(self, table):
        values = table[self.fields]
        colors = self.cmap(self.norm(values), bytes=True)
        colors_hex = rgb_to_hex(colors)
        table[self.color_column] = colors_hex

        return table

    @param.depends('cmap', 'fields', 'norm', watch=True)
    def _updated(self):
        self.updated = True


class PeptideLayoutTransform(Transform):
    """
    Takes data with peptide start, end positions and transforms it to bottom, left, right, top positions of rectangles.

    By default, for a given interval start, end (inclusive, exclusive) a rectangle is returned with a height of spanning from
    (start - 0.5) to (end - 0.5)

    """

    left = param.String('start', doc="Field name to use for left coordinate")
    right = param.String('end', doc="Field name to use for for the right coordinate")

    height = param.Integer(1, doc="Height of the rectangles", constant=True)
    value = param.String('', doc="Optional field name to pass through as value column")
    passthrough = param.List([], doc="Optional field names to pass through (for hovertools)")

    wrap = param.Integer(doc="Amount of peptides to plot on the y axis before wrapping around to y axis top")
    step = param.Integer(5, bounds=(1, None), doc="Step size used for finding 'wrap' when its not specified")
    margin = param.Integer(4, doc="Margin space to keep between peptides when finding 'wrap'")

    transform_type = 'peptide_layout'

    def apply(self, table):
        if not self.wrap:
            self.wrap = autowrap(table[self.left], table[self.right], margin=self.margin, step=self.step)

        # Order of columns determines their role, not the names
        columns = ['x0', 'y0', 'x1', 'y1']  # bottom-left (x0, y0) and top right (x1, y1)
        output_table = pd.DataFrame(index=table.index, columns=columns)
        output_table['x0'] = table[self.left] - 0.5
        output_table['x1'] = table[self.right] - 0.5
        cycle = itertools.cycle(range(self.height*self.wrap, 0, -self.height))  # Starts at y top, cyles with wrap
        yvals = np.array(list(itertools.islice(cycle, len(table)))) # Repeat cycle until the correct length
        output_table['y0'] = yvals - self.height
        output_table['y1'] = yvals

        if self.value:
            output_table['value'] = table[self.value]

        if self.passthrough:
            for item in self.passthrough:
                assert item not in ['value', 'index'], "Invalid field name, 'index' and 'value' names are reserved"
                output_table[item] = table[item]
        output_table['index'] = table.index

        return output_table


class RemoveValueTransform(Transform):
    """Removes entires where values in specified field match specified value"""

    field = param.String(doc='Target field')
    value = param.Number(doc='Rows with this value to remove')

    @property
    def query(self):
        return f'{self.field} != {self.value}'

    def apply(self, table):
        try:
            return table.query(self.query)
        except pd.core.computation.ops.UndefinedVariableError:
            return table