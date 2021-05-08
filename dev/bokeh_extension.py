import numpy as np

from bokeh.core.properties import Float, Bool
from bokeh.models import LayoutDOM
from bokeh.util.compiler import TypeScript
from bokeh.plotting import curdoc
from bokeh.layouts import column
from bokeh.models import Button, CustomJS
import panel as pn
from panel.widgets.base import Widget
import param
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
  class Stage {
    viewerControls: ViewerControls
    constructor(elementId: String, params?: object)
    loadFile(s: String, params?: object):  Promise<StructureComponent>
    autoView() : void
    setSpin(flag: Boolean): void
  }

  class ViewerControls {
        position: Vector3
        Orientation: Matrix4
    }

  class Matrix4{
    elements: Array<Number>
  }

  class Vector3{
    x: number
    y: number
    z: number
  }


  class StructureComponent{
    constructor(stage: Stage, params?: object)
    addRepresentation(type: String, params?: object) : RepresentationComponent
    autoView() : void
  }

  class RepresentationComponent{
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
    console.log(this._stage != null)
    this._stage.loadFile("rcsb://1CRN").then(function (o) {
        o.addRepresentation("cartoon");
        o.autoView();  
        console.log('hi')
    });
    var m = this.model
    var stage = this._stage
    stage.setSpin(m.spin)
    document.addEventListener('click', function(){
       stage.setSpin(m.spin) 
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
    x: p.Property<number>
    y: p.Property<number>
    z: p.Property<number>
    spin: p.Property<boolean>
    representation: p.Property<string>
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
    this.define<ngl.Props>(({Number, Boolean}) => ({
      x:            [ Number ],
      y:            [ Number ],
      z:            [ Number ],
      spin:         [Boolean, false], 
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
    x = Float

    y = Float

    z = Float

    spin = Bool

class NGLview(Widget):
    # Set the Bokeh model to use
    _widget_type = ngl

    # Rename Panel Parameters -> Bokeh Model properties
    # Parameters like title that does not exist on the Bokeh model should be renamed to None
    _rename = {
        "title": None,
    }

    # Parameters to be mapped to Bokeh model properties
    spin = param.Boolean(default=False)
    x = param.Number(0.0, precedence=0)
    y = param.Number(0.0, precedence=0)
    z = param.Number(0.0, precedence=0)


view = ngl(width= 600,height = 600)
but = Button()
but.js_on_click(CustomJS(args=dict(view = view), code = """
    view.spin = !view.spin;
    """))
l = pn.Column(pn.pane.Bokeh(column(view,but)))
l.servable()

