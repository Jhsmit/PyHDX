import param
from panel.widgets.input import Widget, _BkTextInput, StaticText
from panel.widgets import Spinner
from panel.util import as_unicode
from panel.pane import HTML, Markdown
import panel as pn
from dask.distributed import as_completed
from bokeh.models import LayoutDOM
from bokeh.core.properties import String, Bool, List
from bokeh.models import LayoutDOM
from bokeh.util.compiler import TypeScript
import pathlib
from pyhdx.web.bokeh_extensions.ngl_viewer import ngl


class NumericInput(pn.widgets.input.Widget):
    """
    NumericInput allows input of floats with bounds
    """

    type = param.ClassSelector(default=None, class_=(type, tuple),
                               is_instance=True)

    value = param.Number(default=None)

    start = param.Number(default=None, allow_None=True)

    end = param.Number(default=None, allow_None=True)

    _rename = {'name': 'title', 'type': None, 'serializer': None, 'start': None, 'end': None}

    _source_transforms = {'value': """JSON.parse(value.replace(/'/g, '"'))"""}

    _target_transforms = {'value': """JSON.stringify(value).replace(/,/g, ", ").replace(/:/g, ": ")"""}

    _widget_type = _BkTextInput

    def __init__(self, **params):

        super(NumericInput, self).__init__(**params)
        self._state = ''
        self._validate(None)
        self._callbacks.append(self.param.watch(self._validate, 'value'))

    def _validate(self, event):
        if self.type is None: return
        new = self.value
        if not isinstance(new, self.type) and new is not None:
            if event:
                self.value = event.old
            types = repr(self.type) if isinstance(self.type, tuple) else self.type.__name__
            raise ValueError('LiteralInput expected %s type but value %s '
                             'is of type %s.' %
                             (types, new, type(new).__name__))

    def _bound_value(self, value):
        if self.start is not None:
            value = max(value, self.start)
        if self.end is not None:
            value = min(value, self.end)
        return value

    def _process_property_change(self, msg):
        if 'value' in msg and msg['value'] is not None:
            try:
                value = float(msg['value'])
                msg['value'] = self._bound_value(value)
                if msg['value'] != value:
                    self.param.trigger('value')
            except ValueError:
                msg.pop('value')
        if 'placeholder' in msg and msg['placeholder'] is not None:
            try:
                msg['placeholder'] = self._format_value(float(msg['placeholder']))
            except ValueError:
                msg.pop('placeholder')
        return msg

    def _process_param_change(self, msg):
        msg = super(NumericInput, self)._process_param_change(msg)

        if 'start' in msg:
            start = msg.pop('start')
            self.param.value.bounds[0] = start
        if 'end' in msg:
            end = msg.pop('end')
            self.param.value.bounds[1] = end

        if 'value' in msg:
            value = '' if msg['value'] is None else msg['value']
            value = as_unicode(value)
            msg['value'] = value
        msg['title'] = self.name
        return msg


class ColoredStaticText(StaticText):
    _format = '<b>{title}</b>: <span class ="panel-colored-statictext">{value}</span>'


class HTMLTitle(HTML):
    title = param.String(
        doc="""Title"""
    )
    priority = 0
    _rename = dict(HTML._rename, title=None)

    def __init__(self, **params):
        super().__init__(**params)
        self._update_title()

    @param.depends('title', watch=True)
    def _update_title(self):
        self.object = f"""<a class="title" href="" >{self.title}</a>"""


class NGL(Widget):

    _widget_type = ngl

    _rename = {
        "title": None,
    }

    pdb_string = param.String(
        doc="""Raw string of PDB file representing molecular structure to visualize."""
    )
    spin = param.Boolean(
        default=False,
        doc="""Toggle spinning of the molecular structure."""
    )
    representation = param.Selector(
        default='cartoon',
        objects=['ball+stick', 'backbone', 'ball+stick', 'cartoon', 'hyperball', 'licorice',
                 'ribbon', 'rope', 'spacefill', 'surface'],
        doc="""The type of representation used to visualize the molecular structure."""
    )
    rcsb_id = param.String(default="rcsb://1CRN")
    color_list = param.List(default=[["white", "*"]])

    def __init__(self, **params):
        super().__init__(**params)

        self.jscallback(representation="document.dispatchEvent(new Event('representation'));")
        self.jscallback(spin="document.dispatchEvent(new Event('spin'));")
        self.jscallback(rcsb_id="document.dispatchEvent(new Event('rcsb_id'));")
        self.jscallback(color_list="document.dispatchEvent(new Event('color_list'));")
        self.jscallback(pdb_string="document.dispatchEvent(new Event('pdb_string'));")


class LoggingMarkdown(Markdown):
    def __init__(self, header, **params):
        super(LoggingMarkdown, self).__init__(**params)
        self.header = header
        self.contents = ''
        self.object = self.header + self.contents

    def write(self, line):
        self.contents = line + self.contents
        self.object = self.header + self.contents


class ASyncProgressBar(param.Parameterized):
    completed = param.Integer(default=0)
    num_tasks = param.Integer(default=10, bounds=(1, None))

    async def run(self, futures):
        async for task in as_completed(futures):
            with pn.io.unlocked():
                self.completed += 1

    @property
    def value(self):
        value = int(100 * (self.completed / self.num_tasks))
        return max(0, min(value, 100)) # todo check why this is somethings out of bounds

    def reset(self):
        self.completed = 0

    def increment(self):
        self.completed += 1

    @param.depends('completed', 'num_tasks')
    def view(self):
        if self.value:
            return pn.widgets.Progress(active=True, value=self.value, align="center", sizing_mode="stretch_width")
        else:
            return None
