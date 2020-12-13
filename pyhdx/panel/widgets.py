import param
from panel.widgets.input import Widget, _BkTextInput, StaticText
from panel.widgets import Spinner
from panel.util import as_unicode
from panel.pane import HTML, Markdown
import panel as pn
from dask.distributed import as_completed


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


class NGLViewer(HTML):
    pdb_string = param.String(
        doc="""Raw string of PDB file representing molecular structure to visualize."""
    )
    rcsb_id = param.String(
        doc="""ID of PDB structure to fetch from the RCSB Protein Data Bank and visualize."""
    )
    no_coverage = param.Color(
        default='#8c8c8c',
        doc="""Hexadecimal color code to use for residues without coverage"""
    )
    color_list = param.List(default=[], doc="""List of """)
    representation = param.Selector(
        default='cartoon',
        objects=['ball+stick', 'backbone', 'ball+stick', 'cartoon', 'hyperball', 'licorice',
                 'ribbon', 'rope', 'spacefill', 'surface'],
        doc="""The type of representation used to visualize the molecular structure."""
    )
    spin = param.Boolean(
        default=False,
        doc="""Toggle spinning of the molecular structure."""
    )
    priority = 0
    _rename = dict(HTML._rename, pdb_string=None, rcsb_id=None, representation=None, spin=None, color_list=None,
                   no_coverage=None)

    def __init__(self, **params):
        super().__init__(**params)
        self.load_string = \
        f"""
        stage = new NGL.Stage("viewport");
        window.addEventListener( "resize", function( event ){{
            stage.handleResize();
        }}, false );
        stage.loadFile("")""" # this currently gives an error, load empty pdb file?
        self._update_object_from_parameters()

    @property
    def color_array(self):
        """return a string to put into javascript to define the color array"""
        js_string = ', '.join(elem.replace('#', '0x') for elem in self.color_list)
        js_string = js_string.replace('nan', 'noCoverage')
        return js_string

    @param.depends('representation', 'spin', 'color_list', 'no_coverage', watch=True)
    def _update_object_from_parameters(self):
        html =\
            f"""
            <div id="viewport" style="width:100%; height:100%;"></div>
            <script>
            var noCoverage = {self.no_coverage.replace('#', '0x')};
            var colorArray = [{self.color_array}];
            var customScheme = NGL.ColormakerRegistry.addScheme(function (params) {{
                this.atomColor = function (atom) {{
                    if (atom.resno - 1 < colorArray.length) {{
                        return colorArray[atom.resno - 1]
                    }} else {{
                        return noCoverage
                    }}
                }}
            }})

            {self.load_string}.then(function(o){{
                o.addRepresentation("{self.representation}", {{color: customScheme}});
                o.autoView();
                }}
            );
            stage.setSpin({'true' if self.spin else 'false'});
            </script>
            """
        
        self.object = html

    @param.depends('pdb_string', watch=True)
    def _update_object_from_pdb_string(self):
        self.load_string = \
            f"""
            var PDBString = `{self.pdb_string}`;
            stage = new NGL.Stage("viewport");
            window.addEventListener( "resize", function( event ){{
                stage.handleResize();
            }}, false );
            stage.loadFile( new Blob([PDBString], {{type: 'text/plain'}}), {{ ext:'pdb'}} )"""
        self._update_object_from_parameters()

    @param.depends('rcsb_id', watch=True)
    def _update_object_from_rcsb_id(self):
        self.load_string = \
            f"""
            stage = new NGL.Stage("viewport");
            window.addEventListener("resize", function( event ){{
                stage.handleResize();
            }}, false );
            stage.loadFile("rcsb://{self.rcsb_id}")"""
        self._update_object_from_parameters()

    @staticmethod
    def to_hex_string(hex):
        hex.replace('#', '0x')


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
