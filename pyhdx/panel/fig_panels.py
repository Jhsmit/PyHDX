from .base import BokehFigurePanel, FigurePanel, DEFAULT_RENDERERS, DEFAULT_COLORS, MIN_BORDER_LEFT
from .widgets import NGLViewer, LoggingMarkdown
from .log import setup_md_log
from pyhdx.plot import _bokeh_coverage
from bokeh.plotting import figure, curdoc
from bokeh.layouts import column
from bokeh.models import LabelSet, ColumnDataSource, HoverTool, GlyphRenderer, Span, Rect, Range1d
from bokeh.models.markers import Triangle, Circle, Diamond
import panel as pn
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import logging

import param


class CoverageFigure(BokehFigurePanel):
    title = 'Coverage'
    accepted_tags = ['coverage']

    def __init__(self, *args, **params):
        super(CoverageFigure, self).__init__(*args, **params)

    def draw_figure(self):
        fig = figure(title=None, min_border=0, tools='pan,wheel_zoom,box_zoom,save,reset')
        fig.min_border_left = MIN_BORDER_LEFT
        fig.xaxis.axis_label = 'Residue number'

        return fig

    def render_sources(self, src_dict):  #todo , **render_kwargs
        tooltips = [('Pos', '$x{int}'),
                    ('Index', '@index'),
                    ('Start', '@start (@_start)'),
                    ('End', '@end (@_end)'),
                    ('Sequence', '@sequence'),
                    ('Score', '@scores'),
                    ('Uptake', '@uptake (@uptake_corrected / @ex_residues, @maxuptake)')]

        for name, data_source in src_dict.items():
            glyph = Rect(x='x', y='y', width='width', height=1, fill_color='color')
            renderer = self.figure.add_glyph(data_source.source, glyph)
            self.renderers[name] = renderer

            hovertool = HoverTool(renderers=[renderer], tooltips=tooltips)
            self.figure.add_tools(hovertool)


class LinearLogFigure(BokehFigurePanel):
    """base class for bokeh figure which can switch between log and linear axis
        This is a very slim base class (perhaps it should be a mixin?)
    """

    def _redraw_event(self, event):
        """Redraw the figure with new kwargs passed to figure"""
        kwargs = {event.name: event.new.lower()}
        self.redraw(**kwargs)


class ThdFigure(LinearLogFigure):
    """base class for figures extending LinearLogFigure with threshold lines
    parent must have ClassificationControl
    """
    accepted_tags = ['mapping']

    # def __init__(self, parent, *args, **params):
    #     super(ThdFigure, self).__init__(parent, *args, **params)

    def draw_figure(self, **kwargs):
        y_axis_type = kwargs.pop('y_axis_type', 'log')
        x_axis_type = kwargs.pop('x_axis_type', 'linear')
        fig = figure(y_axis_type=y_axis_type, x_axis_type=x_axis_type,
                     tools='pan,wheel_zoom,box_zoom,save,reset', **kwargs)
        fig.min_border_left = MIN_BORDER_LEFT
        fig.xaxis.axis_label = 'Residue number'
        fig.yaxis.axis_label = self.y_label

        # todo refactor controller access
        for _ in range(self.control_panels['ClassificationControl'].param['num_colors'].bounds[1] - 1):
            sp = Span(location=0, dimension='width')
            sp.tags = ['thd']
            sp.visible = False
            fig.add_layout(sp)

        return fig

    def render_sources(self, src_dict, **render_kwargs):
        for name, data_source in src_dict.items():
            kwargs = {**data_source.render_kwargs, **render_kwargs}
            # try:
            #     y_name = render_kwargs.pop('y')
            # except KeyError:
            #     y_name = data_source.render_kwargs.pop('y')

            glyph_func = getattr(self.figure, data_source.renderer)
            renderer = glyph_func(source=data_source.source, name=name,
                                  legend_label=name, **kwargs)  #todo size is being specified at two different places now

            self.renderers[name] = renderer
            hovertool = HoverTool(renderers=[renderer],
                                  tooltips=[('Residue', '@r_number{int}'), (self.y_label, f'@{kwargs["y"]}')],
                                  mode='vline')
            self.figure.add_tools(hovertool)

        if self.renderers:
            self.figure.legend.click_policy = 'hide'

    def _draw_thds(self, *events):
        #todo duplicate code, subclass
        spans = self.figure.select(tags='thd')

        if not self.figure.renderers:
            return

        y_names, source_names = zip(*[(renderer.glyph.y, renderer.data_source.name) for renderer in self.figure.renderers])
        if not self.control_panels['ClassificationControl'].target in source_names:
            self._hide_thds()
            return
        if not self.control_panels['ClassificationControl'].quantity in y_names:
            self._hide_thds()
            return

        spans.sort(key=lambda x: x.id)
        for i, span in enumerate(spans):
            if i < len(self.control_panels['ClassificationControl'].values):
                span.location = self.control_panels['ClassificationControl'].values[i]
                span.visible = self.control_panels['ClassificationControl'].show_thds
            else:
                span.visible = False

    def _hide_thds(self):
        spans = self.figure.select(tags='thd')
        for i, span in enumerate(spans):
            span.visible = False

class RateFigure(ThdFigure):
    title = 'Rates'
    accepted_tags = [('mapping', 'rate')]
    y_label = 'Rate (min⁻¹)'

    def setup_hooks(self):
        super().setup_hooks()
        self.control_panels['ClassificationControl'].param.watch(self._draw_thds, ['values', 'show_thds'])

    def render_sources(self, src_dict, **render_kwargs):
        super().render_sources(src_dict, **render_kwargs)


class PFactFigure(ThdFigure):
    title = 'Protection Factors'
    accepted_tags = [('mapping', 'pfact')]
    y_label = 'Protection factor'

    def setup_hooks(self):
        #todo move to superclass?
        super().setup_hooks()
        self.control_panels['ClassificationControl'].param.watch(self._draw_thds, ['values', 'show_thds'])

    def render_sources(self, src_dict, **render_kwargs):
        super().render_sources(src_dict, y='pfact', **render_kwargs)


class DeltaGFigure(ThdFigure):
    title = 'DeltaG'
    accepted_tags = [('mapping', 'deltaG')]
    y_label = 'DeltaG (J/mol)'

    def setup_hooks(self):
        #todo move to superclass?
        super().setup_hooks()
        self.control_panels['ClassificationControl'].param.watch(self._draw_thds, ['values', 'show_thds'])

    def draw_figure(self, **kwargs):
        return super().draw_figure(y_axis_type='linear')

    def render_sources(self, src_dict, **render_kwargs):
        super().render_sources(src_dict, y='deltaG', **render_kwargs)


class BinaryComparisonFigure(ThdFigure):
    title = 'Binary Comparison'
    accepted_tags = [('comparison', 'mapping')]  # [ ('x' AND 'y') OR 'z' OR 'asdf']
    x_label = 'Residue number'  # move these to _redraw_kwargs?
    y_label = 'Difference'

    def __init__(self, parent, *args, **params):
        super(BinaryComparisonFigure, self).__init__(parent, *args, **params)
        self._redraw_kwargs['y_axis_type'] = 'linear'

    @property
    def y_kwarg(self):
        try:
            plot_quantity = self.control_panels['ClassificationControl'].quantity
            if plot_quantity is not None:
                return plot_quantity
        except KeyError:
            pass
        return None

    def render_sources(self, src_dict, **render_kwargs):
        kwargs = {**self._render_kwargs, **render_kwargs}
        super().render_sources(src_dict, **kwargs)

    def setup_hooks(self):
        super().setup_hooks()
        self.control_panels['ClassificationControl'].param.watch(self._draw_thds, ['values', 'show_thds'])
        self.control_panels['ClassificationControl'].param.watch(self._log_space_updated, ['log_space'])
        #self.control_panels['ClassificationControl'].param.watch(self._target_updated, ['target'])
        self.control_panels['ClassificationControl'].param.watch(self._quantity_updated, ['quantity'])


    #todo group into one function?
    def _quantity_updated(self, event):
        self._render_kwargs['y'] = event.new
        self.y_label = event.new
        self.redraw(**self._redraw_kwargs) # todo auto pass kwargs?

    # def _target_updated(self, event):  #event. cls, name, new, obj, old, type, what
    #     #todo make redraw accept events
    #     self._render_kwargs['target'] = event.new
    #     self.redraw(**self._redraw_kwargs) # todo auto pass kwargs?

    def _log_space_updated(self, event):
        if event.new:   # True, log space
            self._redraw_kwargs['y_axis_type'] = 'log'
        else:
            self._redraw_kwargs['y_axis_type'] = 'linear'
        self.redraw(**self._redraw_kwargs) # todo auto pass kwargs?

    def draw_figure(self, **kwargs):
        y_axis_type = kwargs.pop('y_axis_type', 'linear')
        return super().draw_figure(y_axis_type=y_axis_type, **kwargs)


class ScoresFigure(LinearLogFigure):
    title = 'Scores'
    accepted_tags = [('scores', 'mapping')]
    x_label = 'Residue number'
    y_label = 'Value'

    def render_sources(self, src_dict):
        for name, data_source in src_dict.items():
            for y_field in data_source.render_kwargs['y']:
                glyph_func = getattr(self.figure, data_source.renderer)
                kwargs = data_source.render_kwargs.copy()
                kwargs.pop('y')

                renderer = glyph_func(**kwargs, y=y_field, source=data_source.source, name=name,
                                      legend_label=f'{y_field}')

                self.renderers[name] = renderer
                # hovertool = HoverTool(renderers=[renderer],
                #                       tooltips=[('Residue', '@r_number{int}'), (self.y_label, f'@{data_source.render_kwargs["y"]}')],
                #                       mode='vline')
                # self.figure.add_tools(hovertool)

            if self.renderers:
                self.figure.legend.click_policy = 'hide'


class SingleValueFigure(LinearLogFigure):
    title = 'Values'
    accepted_tags = [('comparison', 'mapping')]
    x_label = 'Residue number'
    y_label = 'Value'

    def render_sources(self, src_dict):
        for name, data_source in src_dict.items():
            for field, render_func in zip(['value1', 'value2'], ['triangle', 'square']):
                glyph_func = getattr(self.figure, render_func)
                kwargs = data_source.render_kwargs.copy()
                kwargs.pop('y')

                renderer = glyph_func(**kwargs, y=field, source=data_source.source, name=name,
                                      legend_label=name + f'_{field}')

                self.renderers[name] = renderer
                hovertool = HoverTool(renderers=[renderer],
                                      tooltips=[('Residue', '@r_number{int}'), (self.y_label, f'@{data_source.render_kwargs["y"]}')],
                                      mode='vline')
                self.figure.add_tools(hovertool)

            if self.renderers:
                self.figure.legend.click_policy = 'hide'

    def draw_figure(self, **kwargs):
        y_axis_type = kwargs.pop('y_axis_type', 'linear')
        return super().draw_figure(y_axis_type=y_axis_type, **kwargs)


class SingleFigure(BinaryComparisonFigure):
    title = 'Values'
    y_label = 'Value'


class FitResultFigure(BokehFigurePanel):
    title = 'Fit Result'
    accepted_tags = ['uptake_curve']
    y_label = 'Uptake corrected'
    x_label = 'Time'

    def __init__(self, parent, *args, **params):
        super(FitResultFigure, self).__init__(parent, *args, **params)

    def setup_hooks(self):
        self.control_panels['FitResultControl'].param.watch(self._redraw_event, ['x_axis_type'])

    def _redraw_event(self, event):
        """Redraw the figure with new kwargs passed to figure"""
        kwargs = {event.name: event.new.lower()}
        self.redraw(**kwargs)

    def draw_figure(self, **kwargs):
        fig = super().draw_figure(x_axis_type=self.control_panels['FitResultControl'].x_axis_type.lower())
        return fig

    def render_sources(self, src_dict, **render_kwargs):
        super().render_sources(src_dict)

        if self.renderers:
            self.figure.legend.location = "bottom_right"


class ProteinFigure(FigurePanel):
    title = 'Protein View'
    accepted_tags = ['mapping']

    js_files = {'ngl': "https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.37/dist/ngl.js"}

    def __init__(self, *args, **params):
        super(ProteinFigure, self).__init__(*args, **params)
        self.ngl_view = NGLViewer(sizing_mode='stretch_both')

    def setup_hooks(self):
        params = ['rcsb_id', 'no_coverage', 'representation', 'spin']
        self.parent.control_panels['ProteinViewControl'].param.watch(self._update_event, params)
        self.parent.control_panels['ProteinViewControl'].file_input.param.watch(self._update_pdb_file, 'value')
        self.parent.control_panels['ProteinViewControl'].param.watch(self._parent_sources_updated, 'target_dataset')

    def _parent_sources_updated(self, *events):
        target = self.parent.control_panels['ProteinViewControl'].target_dataset
        accepted_sources = {k: src for k, src in self.parent.sources.items() if src.resolve_tags(self.accepted_tags)}
        accepted_sources = {k: src for k, src in accepted_sources.items() if k == target}
        new_items = {k: v for k, v in accepted_sources.items() if k not in self.renderers}
        self.add_sources(new_items)

        removed_items = self.renderers.keys() - self.parent.sources.keys()  # Set difference
        #self.parent.logger.debug(f'Protein view removed items: {removed_items}')
        self.remove_sources(removed_items)

    def _data_updated_callback(self, attr, old, new):
        self._update_colors(new['r_number'], new['color'])

    def render_sources(self, src_dict):
        for name, data_source in src_dict.items():
            self._update_colors(data_source.source.data['r_number'], data_source.source.data['color'])

    @property
    def no_coverage(self):
        return self.parent.control_panels['ProteinViewControl'].no_coverage

    def _update_pdb_file(self, event):
        self.ngl_view.pdb_string = event.new.decode()

    def _update_colors(self, r_number, color_arr):
        r_start = r_number[0]
        color_list = list(color_arr)
        if r_start < 1:
            remove_num = 1 - r_start
            color_list = color_list[remove_num:]
        else:
            fill_num = r_start - 1
            color_list = fill_num*[self.parent.control_panels['ProteinViewControl'].no_coverage] + color_list

        self.ngl_view.color_list = color_list

    def _update_event(self, event):
        setattr(self.ngl_view, event.name, event.new)

    @property
    def panel(self):
        return self.ngl_view


class LoggingFigure(FigurePanel):
    title = 'Logging'

    def __init__(self, *args, **params):
        super(LoggingFigure, self).__init__(*args, **params)
        self.markdown = LoggingMarkdown('### Log Window \n', sizing_mode='stretch_both')

        sh = logging.StreamHandler(self.markdown)
        sh.terminator = '  \n'
        formatter = logging.Formatter('%(asctime)s [%(levelname)s]: %(message)s', "%Y-%m-%d %H:%M:%S")
        sh.setFormatter(formatter)
        #sh.setLevel(logging.DEBUG)
        self.parent.logger.addHandler(sh)

    def setup_hooks(self):
        # todo add function on base class for doing these things (with try/except) ?
        # Or the users has to be responsible and set this up correctly
        # if the hook is unwanted, the class should be subclassed with override on setup_hooks

        self.parent.control_panels['OptionsControl'].param.watch(self._update_log_level, ['log_level'])
        self.parent.control_panels['OptionsControl'].param.trigger('log_level')

    def _update_log_level(self, event):
        self.parent.logger.setLevel(event.new)

    @property
    def panel(self):
        return self.markdown
