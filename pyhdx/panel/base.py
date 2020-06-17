import param
import panel as pn
from bokeh.plotting import figure

DEFAULT_RENDERERS = {'fit1': 'triangle', 'fit2': 'circle'}
DEFAULT_COLORS = {'fit1': 'blue', 'fit2': 'red'}
DEFAULT_CLASS_COLORS = ['#cc0c49', '#eded0e', '#1930e0']



class PanelBase(param.Parameterized):
    """base class for mixin panels"""

    position = ''
    @property
    def panel(self):
        return None


class FigurePanel(PanelBase):
    """"base class for figures"""

    _controlled_by = []  # list of panel controllers

    def __init__(self, parent, controllers, **params):
        self.parent = parent  #main controller
        self.controllers = controllers  #side controllers
        super(FigurePanel, self).__init__(**params)

    def draw_figure(self):
        return figure()

    def _update(self):
        """redraw the graph"""
        self.bk_pane.param.trigger('object')

    @property
    def panel(self):
        return self.bk_pane


class ControlPanel(PanelBase):
    """base class for left control pannels"""

    header = 'Default Header'

    def __init__(self, parent, **params):
        self.parent = parent
        super(ControlPanel, self).__init__(**params)

        self._widget_dict = self.make_dict()
        self._widget_list = self.make_list()  # this list after its made isnt / shouldnt be used?
        self._box = self.make_box()

    def make_box(self):
        md = pn.pane.Markdown(f'### {self.header}')
        return pn.WidgetBox(*([md] + self._widget_list))

    def generate_widgets(self, **kwargs):
        """returns a dict with keys parameter names and values default mapped widgets"""
        return {k: v for k, v in zip(list(self.param)[1:], pn.Param(self.param, show_name=False, widgets=kwargs))}

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

#get_widget = pn.Param.get_widget

