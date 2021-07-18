import numpy as np

from collections import OrderedDict
import bokeh
from bokeh.core.properties import Bool, String, List
from bokeh.models import LayoutDOM
from bokeh.util.compiler import TypeScript
from bokeh.plotting import curdoc
from bokeh.layouts import column
from bokeh.models import Button, CustomJS
import panel as pn
from panel.widgets.base import Widget
import param
from pathlib import Path
from pyhdx.web.base import STATIC_DIR


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
    var first_scheme = NGL.ColormakerRegistry.addSelectionScheme(m.color_list, "new scheme");

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
        stage.loadFile(m.rscb).then(function (o) {
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
    rscb: p.Property<string>
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
      rscb:             [String],
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
    rscb = String
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
    rscb = param.String(default="rcsb://1CRN")
    no_coverage = param.String(default='0x8c8c8c')
    color_list = param.List(default=[["red", "64-74 or 134-154 or 222-254 or 310-310 or 322-326"],
                                     ["green", "311-322"],
                                     ["yellow",
                                      "40-63 or 75-95 or 112-133 or 155-173 or 202-221 or 255-277 or 289-309"],
                                     ["blue", "1-39 or 96-112 or 174-201 or 278-288"],
                                     ["white", "*"]])


class test(param.Parameterized):
    pdb_string = param.String()
    rcsb_id = param.String(default="1qyn", doc='RCSB PDB identifier of protein entry to download and visualize.')
    representation = param.Selector(default='cartoon',
                                    objects=['backbone', 'ball+stick', 'cartoon', 'hyperball', 'licorice',
                                             'ribbon', 'rope', 'spacefill', 'surface'],
                                    doc='Representation to use to render the protein.')
    spin = param.Boolean(default=False, doc='Rotate the protein around an axis.')
    #no_coverage = param.Color(default='#8c8c8c', doc='Color to use for regions of no coverage.')
    color_list = param.List(default=[["red", "64-74 or 134-154 or 222-254 or 310-310 or 322-326"],
                                     ["green", "311-322"],
                                     ["yellow",
                                      "40-63 or 75-95 or 112-133 or 155-173 or 202-221 or 255-277 or 289-309"],
                                     ["blue", "1-39 or 96-112 or 174-201 or 278-288"],
                                     ["white", "*"]])


class watcher(object):
    def __init__(self, ngl_viewer, parent):
        self.ngl_viewer = ngl_viewer
        self.parent = parent
        self.parent.param.watch(self.changeSpin, "spin")
        self.parent.param.watch(self.changeRepresentation, "representation")
        self.parent.param.watch(self.changerscb, "rcsb_id")
        #self.parent.param.watch(self.changeColor, "no_coverage")
        self.parent.param.watch(self.changeColorList, "color_list")
        self.parent.param.watch(self.changePDBString, "pdb_string")


    def changeSpin(self, event):
        self.ngl_viewer.spin = event.new


    def changeRepresentation(self, event):
        print(event.new)
        self.ngl_viewer.representation = event.new


    def changerscb(self, event):
        self.ngl_viewer.rscb = event.new


    def changeColor(self, event):
        self.ngl_viewer.color = event.new.no_coverage.replace('#', '0x')


    def changeColorList(self, event):
        self.ngl_viewer.color_list = event.new

    def changePDBString(self, event):
        print(event.new)
        self.ngl_viewer.pdb_string = event.new


class NGLView_factory:

    @staticmethod
    def create_view():
        view = NGLview(sizing_mode='stretch_both', representation='cartoon')
        view.jscallback(representation="document.dispatchEvent(new Event('representation'));")
        view.jscallback(spin="document.dispatchEvent(new Event('spin'));")
        view.jscallback(rscb="document.dispatchEvent(new Event('rcsb_id'));")
        #view.jscallback(no_coverage="document.dispatchEvent(new Event('no_coverage'));")
        view.jscallback(color_list="document.dispatchEvent(new Event('color_list'));")
        view.jscallback(pdb_string="document.dispatchEvent(new Event('pdb_string'));")
        return view


from bokeh.settings import settings
print(settings.minified)

import os
os.environ["BOKEH_XSRF_COOKIES"] = "True"


p = test()
view = NGLView_factory.create_view()
watch = watcher(view, p)

result = pn.Row(p.param, view)

pdb_string = Path('1qyn.pdb').read_text()
result.servable()

def cb():
    view.pdb_string = pdb_string

pn.state.onload(cb)

current_dir = Path(__file__).parent

#pn.serve(result, static_dirs={'pyhdx': STATIC_DIR, 'bk': current_dir})

# if __name__ == '__main__':
#     p = test()
#     view = NGLView_factory.create_view()
#     watch = watcher(view, p)
#
#     result = pn.Column(p.param, view)
#     pn.serve(result)
