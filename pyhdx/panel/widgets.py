import param
from panel.widgets.input import Widget, _BkTextInput, StaticText
from panel.widgets import Spinner
from panel.pane import HTML


class NumericInput(Widget):
    value = param.Number(default=0, allow_None=True, bounds=[None, None])

    placeholder = param.Number(default=None)

    start = param.Number(default=None, allow_None=True)

    end = param.Number(default=None, allow_None=True)

    _widget_type = _BkTextInput

    formatter = param.Parameter(default=None)

    _rename = {'name': 'title', 'formatter': None, 'start': None, 'end': None}

    def __init__(self, **params):
        if params.get('value') is None:
            value = params.get('start', self.value)
            if value is not None:
                params['value'] = value
        super(NumericInput, self).__init__(**params)

    def _bound_value(self, value):
        if self.start is not None:
            value = max(value, self.start)
        if self.end is not None:
            value = min(value, self.end)
        return value

    def _format_value(self, value):
        if self.formatter is not None:
            value = self.formatter.format(value)
        else:
            value = str(value)
        return value

    def _process_param_change(self, msg):
        msg.pop('formatter', None)

        if 'start' in msg:
            start = msg.pop('start')
            self.param.value.bounds[0] = start
        if 'end' in msg:
            end = msg.pop('end')
            self.param.value.bounds[1] = end

        if 'value' in msg and msg['value'] is not None:
            msg['value'] = self._format_value(self.value)
        if 'placeholder' in msg and msg['placeholder'] is not None:
            msg['placeholder'] = self._format_value(self.placeholder)
        return msg

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


class ColoredStaticText(StaticText):
    _format = '<b>{title}</b>: <span class ="panel-colored-statictext">{value}</span>'


class NGLViewer(HTML):
    pdb_string = param.String()
    rcsb_id = param.String()
    no_coverage = param.Color(default='#8c8c8c')
    color_list = param.List([])
    representation = param.Selector(default='cartoon',
                                    objects=['ball+stick', 'backbone', 'ball+stick', 'cartoon', 'hyperball', 'licorice',
                                             'ribbon', 'rope', 'spacefill', 'surface'])
    spin = param.Boolean(default=False)
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
        stage.loadFile("")""" # this currently gives an error
        self._update_object_from_parameters()

    @param.depends('representation', 'spin', 'color_list', 'no_coverage', watch=True)
    def _update_object_from_parameters(self):
        html =\
            f"""
            <div id="viewport" style="width:100%; height:100%;"></div>
            <script>
            var noCoverage = {self.no_coverage.replace('#', '0x')};
            var colorArray = [{', '.join(elem.replace('#', '0x') for elem in self.color_list)}];
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
        # print('000000000000000')
        # print(html)
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