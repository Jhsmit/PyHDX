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
from bokeh.models import CustomJS


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

TS_CODE = """
// This custom model wraps one part of the third-party vis.js library:
//
//     http://visjs.org/index.html
//
// Making it easy to hook up python data analytics tools (NumPy, SciPy,
// Pandas, etc.) to web presentations using the Bokeh server.

import {LayoutDOM, LayoutDOMView} from "models/layouts/layout_dom"
import {LayoutItem} from "core/layout"
import * as p from "core/properties"

declare namespace NGL {
  class AtomProxy{

  }
  class Blob{
    constructor(list: Array<String>, ob: object)
  }
  class Colormaker{
    atomColor: (atom: AtomProxy) => string
  }

  class ColormakerRegistry{
    static addScheme(scheme: () => void) : String
    static addSelectionScheme(dataList: Array<Array<String>>, label: String): String
  }

  class Component{
    removeAllRepresentations(): void
    addRepresentation(type: String, params?: object) : RepresentationComponent
    reprList: RepresentationElement[]
  }

  class Matrix4{
    elements: Array<Number>
  }

  class RepresentationComponent{
  }

  class RepresentationElement{
    setParameters(params: any): this
    getParameters(): object
  }

  class Stage {
    compList: Array<Component>
    viewerControls: ViewerControls
    constructor(elementId: String, params?: object)
    loadFile(s: String| Blob, params?: object):  Promise<StructureComponent>
    autoView() : void
    setSpin(flag: Boolean): void
    removeAllComponents(type: String): void
    addRepresentation(representation: String): void
  }

  class ScriptComponent{
    constructor(stage: Stage, params?: object)
    addRepresentation(type: String, params?: object) : RepresentationComponent
    autoView() : void
    removeAllRepresentations(): void
    reprList: RepresentationElement[]
  }

  class StructureComponent{
    constructor(stage: Stage, params?: object)
    addRepresentation(type: String, params?: object) : RepresentationComponent
    autoView() : void
    removeAllRepresentations(): void
    reprList: RepresentationElement[]
  }

  class SurfaceComponent{
    constructor(stage: Stage, params?: object)
    addRepresentation(type: String, params?: object) : RepresentationComponent
    autoView() : void
    removeAllRepresentations(): void
    reprList: RepresentationElement[]
  }



  class Vector3{
    x: number
    y: number
    z: number
  }

  class ViewerControls {
        position: Vector3
        Orientation: Matrix4
    }
}

export class NGLView extends LayoutDOMView {
  model: ngl
  public spin: Boolean
  public _stage: NGL.Stage

  initialize(): void {
    super.initialize()

    const url = "https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.37/dist/ngl.js"
    const script = document.createElement("script")
    script.onload = () => this._init()
    script.async = false
    script.src = url
    document.head.appendChild(script)
  }

  public set_variable_x(x: number): void {
    this._stage.viewerControls.position.x = x;
  }

  private _init(): void {
    // Create a new Graph3s using the vis.js API. This assumes the vis.js has
    // already been loaded (e.g. in a custom app template). In the future Bokeh
    // models will be able to specify and load external scripts automatically.
    //
    // BokehJS Views create <div> elements by default, accessible as this.el.
    // Many Bokeh views ignore this default <div>, and instead do things like
    // draw to the HTML canvas. In this case though, we use the <div> to attach
    // a Graph3d to the DOM.


    this.el.setAttribute('id','viewport')
    console.log("the id is: " + this.el.getAttribute('id'))
    this._stage = new NGL.Stage('viewport')
    var m = this.model
    var stage = this._stage
    var promise = this._stage.loadFile("rcsb://1CRN")
    var first_scheme = NGL.ColormakerRegistry.addSelectionScheme(m.color_list, "new scheme");
    promise.then(function (o) {
        o.addRepresentation("cartoon", { color: first_scheme });
        o.autoView();
    });

    stage.setSpin(m.spin)
    document.addEventListener('spin', function(){
       stage.setSpin(m.spin);
    });

    document.addEventListener('representation', function(){
        stage.compList[0].removeAllRepresentations();
        stage.compList[0].addRepresentation(m.representation, { color: first_scheme })
    });

    document.addEventListener('rcsb_id', function(){
        stage.removeAllComponents("");
        stage.loadFile(m.rcsb_id).then(function (o) {
            o.addRepresentation(m.representation, { color: first_scheme })
            o.autoView()
        });
    });

    document.addEventListener('color_list', function(){
        console.log(m.color_list)
        var list: Array<Array<String>> = m.color_list
        try{
              var new_scheme = NGL.ColormakerRegistry.addSelectionScheme(list, "new scheme");
              stage.compList[0].reprList[0].setParameters( { color: new_scheme } );
        }
        catch(err) {
            console.log("badly defined color")
        }
    });

    document.addEventListener('pdb_string', function(){
        stage.removeAllComponents("");
        stage.loadFile( new Blob([m.pdb_string], {type: 'text/plain'}), { ext:'pdb'}).then(function (o) {
            o.addRepresentation(m.representation, { color: first_scheme })
            o.autoView()
        });

    });
   }

  // This is the callback executed when the Bokeh data has an change. Its basic
  // function is to adapt the Bokeh data source to the vis.js DataSet format.
  //get_data(): vis.DataSet {
  //  const data = new vis.DataSet()
  //  const source = this.model.data_source
  //  for (let i = 0; i < source.get_length()!; i++) {
  //    data.add({
  //      x: source.data[this.model.x][i],
  //      y: source.data[this.model.y][i],
  //      z: source.data[this.model.z][i],
  //    })
  //  }
  //  return data
  //}

  get child_models(): LayoutDOM[] {
    return []
  }

  _update_layout(): void {
    this.layout = new LayoutItem()
    this.layout.set_sizing(this.box_sizing())
  }
}

// We must also create a corresponding JavaScript BokehJS model subclass to
// correspond to the python Bokeh model subclass. In this case, since we want
// an element that can position itself in the DOM according to a Bokeh layout,
// we subclass from ``LayoutDOM``
export namespace ngl {
  export type Attrs = p.AttrsOf<Props>

  export type Props = LayoutDOM.Props & {
    spin: p.Property<boolean>
    representation: p.Property<string>
    rcsb_id: p.Property<string>
    no_coverage: p.Property<string>
    color_list: p.Property<any>
    pdb_string: p.Property<string>
  }
}

export interface ngl extends ngl.Attrs {}

export class ngl extends LayoutDOM {
  properties: ngl.Props
  __view_type__: NGLView

  constructor(attrs?: Partial<ngl.Attrs>){
    super(attrs)
  }

  // The ``__name__`` class attribute should generally match exactly the name
  // of the corresponding Python class. Note that if using TypeScript, this
  // will be automatically filled in during compilation, so except in some
  // special cases, this shouldn't be generally included manually, to avoid
  // typos, which would prohibit serialization/deserialization of this model.
  static __name__ = "ngl"

  static init_ngl() {
    // This is usually boilerplate. In some cases there may not be a view.
    this.prototype.default_view = NGLView

    // The @define block adds corresponding "properties" to the JS model. These
    // should basically line up 1-1 with the Python model class. Most property
    // types have counterparts, e.g. ``bokeh.core.properties.String`` will be
    // ``String`` in the JS implementation. Where the JS type system is not yet
    // as rich, you can use ``p.Any`` as a "wildcard" property type.
    this.define<ngl.Props>(({String, Boolean, Any}) => ({
      spin:             [Boolean, false],
      representation:   [String],
      rcsb_id:          [String],
      no_coverage:      [String],
      color_list:       [Any],
      pdb_string:       [String]
    })
    )
  }
}
"""


# This custom extension model will have a DOM view that should layout-able in
# Bokeh layouts, so use ``LayoutDOM`` as the base class. If you wanted to create
# a custom tool, you could inherit from ``Tool``, or from ``Glyph`` if you
# wanted to create a custom glyph, etc.
class ngl(LayoutDOM):
    # The special class attribute ``__implementation__`` should contain a string
    # of JavaScript code that implements the browser side of the extension model.

    # Below are all the "properties" for this model. Bokeh properties are
    # class attributes that define the fields (and their types) that can be
    # communicated automatically between Python and the browser. Properties
    # also support type validation. More information about properties in
    # can be found here:__implementation__ = TypeScript(TS_CODE)
    #
    #    https://docs.bokeh.org/en/latest/docs/reference/core/properties.html#bokeh-core-properties
    __implementation__ = TypeScript(TS_CODE)
    # This is a Bokeh ColumnDataSource that can be updated in the Bokeh
    # server by Python code

    # The vis.js library that we are wrapping expects data for x, y, and z.
    # The data will actually be stored in the ColumnDataSource, but these
    # properties let us specify the *name* of the column that should be
    # used for each field.

    spin = Bool
    representation = String
    rcsb_id = String
    no_coverage = String
    color_list = List(List(String))
    pdb_string = String


class NGLview(Widget):
    # Set the Bokeh model to use
    _widget_type = ngl

    # Rename Panel Parameters -> Bokeh Model properties
    # Parameters like title that does not exist on the Bokeh model should be renamed to None
    _rename = {
        "title": None,
    }
    pdb_string = param.String()
    # Parameters to be mapped to Bokeh model properties
    spin = param.Boolean(default=False)
    representation = param.String(default="cartoon")
    rcsb_id = param.String(default="rcsb://1CRN")
    no_coverage = param.String(default='0x8c8c8c')
    color_list = param.List(default=[["red", "64-74 or 134-154 or 222-254 or 310-310 or 322-326"],
                                     ["green", "311-322"],
                                     ["yellow",
                                      "40-63 or 75-95 or 112-133 or 155-173 or 202-221 or 255-277 or 289-309"],
                                     ["blue", "1-39 or 96-112 or 174-201 or 278-288"],
                                     ["white", "*"]])


class NGLView_factory:

    @staticmethod
    def create_view(**kwargs):
        view = NGLview(**kwargs)
        view.jscallback(representation="document.dispatchEvent(new Event('representation'));")
        view.jscallback(spin="document.dispatchEvent(new Event('spin'));")
        view.jscallback(rcsb_id="document.dispatchEvent(new Event('rcsb_id'));")
        # view.jscallback(no_coverage="document.dispatchEvent(new Event('no_coverage'));")
        view.jscallback(color_list="document.dispatchEvent(new Event('color_list'));")
        view.jscallback(pdb_string="document.dispatchEvent(new Event('pdb_string'));")
        return view

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
