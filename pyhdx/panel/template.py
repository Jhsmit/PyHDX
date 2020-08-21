import pathlib

from panel.template import GoldenTemplate
import panel as pn


class ExtendedGoldenTemplate(GoldenTemplate):

    _template = pathlib.Path(__file__).parent / 'golden.html'

#    _css = pathlib.Path(__file__).parent / 'golden.css'


class SubString(str):
    """
    Extends the `string` class such that it can be used to monkey-patch the _template class attribute of GoldenTemplate
    """
    def read_into(self):
        return str(self)


class GoldenElvis(object):
    """
    Adaptation of Leon van Kouwen's elvis layout system
    https://github.com/LeonvanKouwen/elvis

    Generates a jinja GoldenLayout Template based on panel's default GoldenLayout Template. This modification features a
    fixed sidebar with a main layout part which can be customized with columns/rows/stacks

    """

    NESTABLE = \
        """
        {
            type: '%s',
            content: [ %s ]
        },
        """

    VIEW = \
        """
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

    def __init__(self, template, theme, title=None):
        
        self.template = template(title=title, theme=theme)

    @property
    def jinja_base(self):
        _base = pathlib.Path(__file__).parent / 'jinja_base.html'

        return _base.read_text()

    def make_sidebar(self, controllers):
        controls = pn.Column(*[controller.panel for controller in controllers])
        self.template.sidebar.append(controls)


    def compose(self, golden_layout_string):
        """
        Creates a servable template from a golden layout js code string.
        :param golden_layout_string: Result of nesting stacks, columns, rows, and panels
                                     using the methods in this class.
        """
        template = self.jinja_base % golden_layout_string
        self.app = pn.Template(template=template)
        for panel_ID, panel in self.panels.items():
            self.app.add_panel(panel_ID, panel)

    def view(self, view, title=None, width=None, height=None, scrollable=True):
        """
        Adds a viewable panel.
        :param view: The panel to show in this golden layout sub section.
        :param title: The text to show at the top of the panel.
        :param width: Initial width.
        :param height: Initial height.
        """

        # We need to register every panel with a unique name such that after
        # composing the jinja2 template, we can add them (see compose function).
        self.counter = self.counter + 1
        panel_ID = "panel_" + str(self.counter)
        self.panels[panel_ID] = pn.panel(view, sizing_mode='stretch_both')
        title_str = "title: '%s'," % str(title) if title is not None else "title: '',"
        width_str = "width: %s," % str(width) if width is not None else ""
        height_str = "height: %s," % str(height) if height is not None else ""
        scroll_str = "css_classes: ['not_scrollable']" if not scrollable else ""
        settings = title_str + height_str + width_str + scroll_str
        return ClientSideCodeStrings.VIEW % (panel_ID, settings)

    def _block(self, *args, type='stack'):
        """
        Creates nestable js code strings. Note that 'stack', 'colum' and 'row' are the
        strings dictated by the golden layout js code.
        """
        content = ''.join(arg for arg in args)
        return self.NESTABLE % (type, content)

    def stack(self, *args):
        """ Adds a 'tab' element."""
        return self._block(*args, type='stack')

    def column(self, *args):
        """ Vertically aligned panels"""
        return self._block(*args, type='column')

    def row(self, *args):
        """ Horizontally aligned panels"""
        return self._block(*args, type='row')