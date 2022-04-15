import os
import pathlib
import string

import panel as pn

from panel.template import GoldenTemplate

from panel.util import url_path
from param.parameterized import default_label_formatter

from pyhdx.web.widgets import HTMLTitle

dist_path = "/pyhdx/"

SIDEBAR_WIDTH = 300


class ExtendedGoldenTemplate(GoldenTemplate):

    pass


class ReadString(str):
    """
    Extends the `string` class such that it can be used to monkey-patch the _template class attribute of GoldenTemplate
    """

    def read_text(self):
        return str(self)


class GoldenElvis(object):
    """
    Adaptation of Leon van Kouwen's elvis layout system
    https://github.com/LeonvanKouwen/elvis

    Generates a jinja GoldenLayout Template based on panel's default GoldenLayout Template. This modification features a
    fixed sidebar with a main layout part which can be customized with columns/rows/stacks

    """

    NESTABLE = """
        {
            type: '%s',
            content: [ %s ],
            %s
        },
        """

    VIEW = """
        {   
            type: 'component',
            componentName: 'view',
            componentState: 
            { 
                model: '{{ embed(roots.%s) }}',
                %s
            },
            isClosable: false,
        },
        """

    def __init__(self, main_controller, template, theme, title=None):
        self.main_controller = main_controller
        self.template_cls = template
        self.theme_cls = theme
        self.title = title

        self.panels = {}

    @property
    def jinja_base_string_template(self):
        _base = pathlib.Path(__file__).parent / "jinja_base.html"
        base_string_template = string.Template(_base.read_text())

        return base_string_template

    def compose(self, golden_layout_string, **kwargs):
        """
        Creates a servable template from a golden layout js code string.
        :param main_controller: Application main controller
        :param golden_layout_string: Result of nesting stacks, columns, rows, and panels
                                     using the methods in this class.
        """

        controllers = self.main_controller.control_panels.values()
        template_code = ReadString(
            self.jinja_base_string_template.substitute(main_body=golden_layout_string)
        )
        self.template_cls._template = template_code

        template = self.template_cls(title=self.title, theme=self.theme_cls, **kwargs)
        controls = pn.Accordion(
            *[controller.panel for controller in controllers],
            toggle=True,
            sizing_mode="fixed",
            width=SIDEBAR_WIDTH,
        )

        template.sidebar.append(controls)

        for panel_ID, panel in self.panels.items():
            template._render_items[panel_ID] = (panel, ["main"])

        return template

    def view(self, view_name, title=None, width=None, height=None, scrollable=True):
        """
        Adds a viewable panel.
        :param view: The panel to show in this golden layout sub section.
        :param title: The text to show at the top of the panel.
        :param width: Initial width.
        :param height: Initial height.
        """
        # pn.config.js_files.update(fig_panel.js_files)

        # We need to register every panel with a unique name such that after
        # composing the jinja2 template, we can add them (see compose function).

        # It seems that these unique names cannot start with a number or they cannot be referenced directly
        # Therefore, currently tmpl.main.append cannot be used as this generates
        fig_panel = self.main_controller.views[view_name]
        panel_ID = "ID" + str(id(fig_panel))
        title = default_label_formatter(title or getattr(fig_panel, "name", None))

        fig_panel.update()  # intialize
        item = pn.Row(
            fig_panel.panel, sizing_mode="stretch_both"
        )  # Place figure in layout
        self.panels[panel_ID] = item
        title_str = "title: '%s'," % str(title) if title is not None else "title: '',"
        width_str = "width: %s," % str(width) if width is not None else ""
        height_str = "height: %s," % str(height) if height is not None else ""
        scroll_str = "css_classes: ['lm_content_noscroll']" if not scrollable else ""

        # scroll_str = "css_classes: ['overflow-y: hidden !important']" # this doesnt work
        # scroll_str = "overflow: 'hidden'," #if not scrollable else ""
        settings = title_str + height_str + width_str + scroll_str
        return self.VIEW % (panel_ID, settings)

    def get_settings(self, **kwargs):
        settings = ""
        for name, val in kwargs.items():
            if isinstance(val, str):
                settings += f"{name}: '{val}'"
            elif isinstance(val, (float, int)):
                settings += f"{name}: {val}"

        return settings

    def _block(self, *args, container="stack", **kwargs):
        """
        Creates nestable js code strings. Note that 'stack', 'colum' and 'row' are the
        strings dictated by the golden layout js code.
        """
        content = "".join(arg for arg in args)
        settings = self.get_settings(**kwargs)
        return self.NESTABLE % (container, content, settings)

    def stack(self, *args, **kwargs):
        """Adds a 'tab' element."""
        return self._block(*args, container="stack", **kwargs)

    def column(self, *args, **kwargs):
        """Vertically aligned panels"""
        return self._block(*args, container="column", **kwargs)

    def row(self, *args, **kwargs):
        """Horizontally aligned panels"""
        return self._block(*args, container="row", **kwargs)
