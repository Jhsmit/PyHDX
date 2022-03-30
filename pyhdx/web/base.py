from pathlib import Path

import panel as pn
import param

from pyhdx.web.main_controllers import MainController
from pyhdx.web.template import SIDEBAR_WIDTH
from pyhdx.web.utils import get_view

DEFAULT_RENDERERS = {
    "half-life": "hex",
    "fit1": "triangle",
    "fit2": "circle",
    "TF_rate": "diamond",
    "pfact": "circle",
}
DEFAULT_COLORS = {
    "half-life": "#f37b21",
    "fit1": "#2926e0",
    "fit2": "#f20004",
    "TF_rate": "#03ab1d",
    "pfact": "#16187d",
    "uptake_corrected": "#000000",
    "fr_pfact": "#ba0912",
}
# DEFAULT_CLASS_COLORS = ['#0e1875', '#fdaf61', '#d73027']  # rigid to flexible
DEFAULT_CLASS_COLORS = ["#0a0ac2", "#0ac20a", "#c20a0a"][::-1]  #  (HSL xxx, 90, 40)
# DEFAULT_CLASS_COLORS = ['#3d3df5', '#3df53d', '#f53d3d'][::-1] #  (HSL xxx, 90, 60)

MIN_BORDER_LEFT = 65
STATIC_DIR = Path(__file__).parent / "static"


class ControlPanel(param.Parameterized):
    """base class for left control pannels"""

    _type = None

    header = "Default Header"

    parent = param.ClassSelector(MainController, precedence=-1)

    _excluded = param.List(
        [],
        precedence=-1,
        doc="parameters whose widgets are excluded from the control panel view. This list can be modified to update "
        "widget layout",
    )

    def __init__(self, parent, **params):
        # self.parent = parent
        # self.ex
        super(ControlPanel, self).__init__(parent=parent, **params)

        self.widgets = (
            self.make_dict()
        )  # atm on some objects this is a list, others dict
        self._box = self.make_box()  # _panel equivalent

        if self._layout:
            for widget_source, contents in self._layout:
                if widget_source != "self":
                    _type, name = widget_source.split(".")
                    if _type == "transforms":
                        object = getattr(self, _type)[name]
                        object.param.watch(self.update_box, ["redrawn"])

    @property  # todo base class
    def own_widget_names(self):
        return [name for name in self.widgets.keys() if name not in self._excluded]

    @property
    def _layout(self):
        return [
            ("self", self.own_widget_names),
        ]

    @property
    def sources(self):
        return self.parent.sources

    @property
    def transforms(self):
        return self.parent.transforms

    @property
    def opts(self):
        return self.parent.opts

    @property
    def views(self):
        return self.parent.views

    def make_box(self):
        return pn.Column(*self.widget_list, name=self.header, width=SIDEBAR_WIDTH)

    def _update_box(self, *events):
        self.update_box()

    def update_box(self, *events):
        self._box[:] = self.widget_list

    def generate_widgets(self, **kwargs):
        """returns a dict with keys parameter names and values default mapped widgets"""

        # todo respect precedence
        names = [
            p
            for p in self.param
            if self.param[p].precedence is None or self.param[p].precedence > 1
        ]
        widgets = pn.Param(
            self.param, show_name=False, show_labels=True, widgets=kwargs
        )

        return {k: v for k, v in zip(names[1:], widgets)}

    @property
    def widget_list(self):
        """

        Example _layout definitions

        Returns
        -------

        """

        try:
            self._layout
        except AttributeError:
            return [get_view(widget) for widget in self.widgets.values()]

        if self._layout is None:
            return [get_view(widget) for widget in self.widgets.values()]
        else:
            widget_list = []
            for widget_source, contents in self._layout:
                if widget_source == "self":
                    object = self
                else:
                    _type, name = widget_source.split(".")
                    object = getattr(self, _type)[name]

                if isinstance(contents, list):
                    for item in contents:
                        widget_list.append(get_view(object.widgets[item]))
                elif isinstance(contents, str):
                    widget_list.append(get_view(object.widgets[contents]))
                elif contents is None:
                    if hasattr(object, "widgets"):
                        for item in object.widgets.values():
                            widget_list.append(get_view(item))
                    else:
                        panel = object.panel
                        if isinstance(panel, pn.layout.ListLike):
                            for item in panel:
                                widget_list.append(get_view(item))
                        else:
                            widget_list.append(get_view(panel))

        return widget_list

    def make_dict(self):
        """dict of widgets to be shown
        override this method to get custom mapping

        """
        return self.generate_widgets()

    def box_index(self, p_name_or_widget):
        "" "return the index of the widget in the box with parameter p_name"
        if isinstance(p_name_or_widget, str):
            return list(self._box).index(self.widget_dict[p_name_or_widget])
        else:
            return list(self._box).index(p_name_or_widget)

    def box_pop(self, p_name_or_widget):
        """remove the widget with parameter name name from the box"""
        index = self.box_index(p_name_or_widget)
        self._box.pop(index)

    def box_insert_after(self, name_or_widget_after, name_or_widget_insert):
        """insert widget corresponding to parameter with name after the widget name_after"""
        index = self.box_index(name_or_widget_after)
        if isinstance(name_or_widget_insert, str):
            widget = self.widget_dict[name_or_widget_insert]
        else:
            widget = name_or_widget_insert
        self._box.insert(index + 1, widget)

    def get_widget(self, param_name, widget_type, **kwargs):
        """get a single widget with for parameter param_name with type widget_type"""

        # not sure if this function still exists
        return pn.Param.get_widget(
            getattr(self.param, param_name), widget_type, **kwargs
        )[0]

    @property
    def panel(self):
        return self._box
