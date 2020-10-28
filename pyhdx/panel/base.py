import param
import panel as pn
from bokeh.plotting import figure
from functools import partial
from pathlib import Path

DEFAULT_RENDERERS = {'half-life': 'hex', 'fit1': 'triangle', 'fit2': 'circle', 'TF_rate': 'diamond', 'pfact': 'circle'}
DEFAULT_COLORS = {'half-life': '#f37b21', 'fit1': '#2926e0', 'fit2': '#f20004', 'TF_rate': '#03ab1d', 'pfact': '#16187d',
                  'uptake_corrected': '#000000', 'fr_pfact': '#ba0912'}
#DEFAULT_CLASS_COLORS = ['#0e1875', '#fdaf61', '#d73027']  # rigid to flexible
DEFAULT_CLASS_COLORS = ['#0a0ac2', '#0ac20a', '#c20a0a']  # rigid to flexible  (HSL xxx, 90, 40)


MIN_BORDER_LEFT = 65
STATIC_DIR = Path(__file__).parent / 'static'


class PanelBase(param.Parameterized):
    """base class for all panel classes"""

    title = ''

    @property
    def panel(self):
        return None


class FigurePanel(PanelBase):
    """Base class for figure panels"""
    accepted_tags = []
    js_files = {}

    def __init__(self, parent, sources=None, **params):
        super(FigurePanel, self).__init__(**params)
        self.parent = parent  # main controller
        self.parent.param.watch(self._parent_sources_updated, ['sources'])

        sources = sources if sources is not None else {}
        self._render_kwargs = {}
        self._redraw_kwargs = {}

        self.renderers = {}  # dict of renderers
        self.data_sources = {}  # dict with DataSource objects
        self.add_sources(sources)

        self.setup_hooks()

    def setup_hooks(self):
        """override to add watchers to controllers"""
        pass

    @property
    def control_panels(self):
        return self.parent.control_panels

    def _parent_sources_updated(self, *events):
        accepted_sources = {k: src for k, src in self.parent.sources.items() if src.resolve_tags(self.accepted_tags)}
        new_items = {k: v for k, v in accepted_sources.items() if k not in self.renderers}
        self.add_sources(new_items)

        removed_items = self.renderers.keys() - self.parent.sources.keys()  # Set difference
        self.remove_sources(removed_items)

    def add_sources(self, src_dict):
        """add a DataSource object to the figure
        """
        #todo check for already in self.data_sources
        for name, data_source in src_dict.items():
            data_source.source.on_change('data', self._data_updated_callback)
            self.data_sources[name] = data_source

        self.render_sources(src_dict)

    def remove_sources(self, names):
        """remove source from renderers dict and figure"""
        #todo not really sure if this works
        for name in names:
            renderer = self.renderers[name]
            renderer.data_source.remove_on_change('data', self._data_updated_callback)
            self.renderers.pop(name)

    def render_sources(self, src_dict):
        """override to customize how sources are rendered"""
        pass

    def _data_updated_callback(self, attr, old, new):
        """overload for custom callback when data updates"""
        pass

    def update(self):
        """called to update the representation"""
        pass

    @property
    def panel(self):
        raise NotImplementedError('panel property not implemented')


class BokehFigurePanel(FigurePanel):
    """
    Base class of bokeh-based figure panels
    """
    x_label = ''
    y_label = ''

    def __init__(self, parent, sources=None, **params):
        super(BokehFigurePanel, self).__init__(parent, sources=sources, **params)
        self.figure = self.draw_figure()
        self.bk_pane = pn.pane.Bokeh(self.figure, sizing_mode='stretch_both', name=self.title)



    def draw_figure(self, **kwargs):
        """Overload to create a custom figure"""

        fig = figure(**kwargs)
        fig.xaxis.axis_label = self.x_label
        fig.yaxis.axis_label = self.y_label

        return fig

    def redraw(self, **kwargs):
        """calls draw_figure to make a new figure and then redraws all renderers"""

        #src_dict = self.data_sources
        #self.remove_sources(src_dict.keys())
        self.renderers = {}
        #self.renderers = {}
        self.figure = self.draw_figure(**kwargs)
        #self.add_sources(src_dict)
        # todo does the old figure linger on?
        self.render_sources(self.data_sources)
        self.bk_pane.object = self.figure

    def _data_updated_callback(self, attr, old, new):
        if self.parent.doc is not None:
            callback = partial(self.bk_pane.param.trigger, 'object')
            self.parent.doc.add_next_tick_callback(callback)
        else:
            self.bk_pane.param.trigger('object')
            
    def render_sources(self, src_dict, **render_kwargs):
        """override to customize how sources are rendered"""
        for name, data_source in src_dict.items():
            glyph_func = getattr(self.figure, data_source.renderer)
            renderer = glyph_func(**data_source.render_kwargs, source=data_source.source, legend_label=name, name=name,
                                  **render_kwargs)
            self.renderers[name] = renderer

    def remove_sources(self, names):
        """remove source from renderers dict and figure"""
        #todo not really sure if this works
        for name in names:
            renderer = self.renderers[name]
            renderer.data_source.remove_on_change('data', self._data_updated_callback)
            self.figure.renderers.remove(renderer)
            for tool in self.figure.tools:
                try:
                    tool.renderers.remove(renderer)
                except (ValueError, AttributeError):
                    pass
            self.renderers.pop(name)

    def update(self):
        callback = partial(self.bk_pane.param.trigger, 'object')
        self.parent.doc.add_next_tick_callback(callback)

    @property
    def panel(self):
        return self.bk_pane


class ControlPanel(PanelBase):
    """base class for left control pannels"""

    header = 'Default Header'

    def __init__(self, parent, **params):
        self.parent = parent
        super(ControlPanel, self).__init__(**params)

        self._widget_dict = self.make_dict()  # should maybe not be private
        self._widget_list = self.make_list()  # this list after its made isnt / shouldnt be used?
        self._box = self.make_box()

    def make_box(self):
        return pn.Card(title=self.header, collapsed=True, *self._widget_list)

    def generate_widgets(self, **kwargs):
        """returns a dict with keys parameter names and values default mapped widgets"""

        names = [p for p in self.param if self.param[p].precedence is None or self.param[p].precedence > 1]
        widgets = pn.Param(self.param, show_name=False, show_labels=True, widgets=kwargs)

        return {k: v for k, v in zip(names[1:], widgets)}
        #
        # return {k: v for k, v in zip(list(self.param)[1:],
        #                              pn.Param(self.param, show_name=False, show_labels=True, widgets=kwargs))}

    def make_list(self):
        """override this method to modify mapping of dict to list"""
        return list(self._widget_dict.values())

    def make_dict(self):
        """dict of widgets to be shown
        override this method to get custom mapping

        """
        return self.generate_widgets()

    def box_index(self, p_name_or_widget):
        ""'return the index of the widget in the box with parameter p_name'
        if isinstance(p_name_or_widget, str):
            return list(self._box).index(self._widget_dict[p_name_or_widget])
        else:
            return list(self._box).index(p_name_or_widget)

    def box_pop(self, p_name_or_widget):
        """remove the widget with parameter name name from the box"""
        index = self.box_index(p_name_or_widget)
        self._box.pop(index)

    def box_insert_after(self, name_or_widget_after, name_or_widget_insert):
        """insert widget corresponding to parameter with name after the widget name_after """
        index = self.box_index(name_or_widget_after)
        if isinstance(name_or_widget_insert, str):
            widget = self._widget_dict[name_or_widget_insert]
        else:
            widget = name_or_widget_insert
        self._box.insert(index + 1, widget)

    def get_widget(self, param_name, widget_type, **kwargs):
        """get a single widget with for parameter param_name with type widget_type"""
        return pn.Param.get_widget(getattr(self.param, param_name), widget_type, **kwargs)[0]

    @property
    def panel(self):
        return self._box


