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
    handleResize(): void
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

    stage.loadFile( new Blob([m.pdb_string], {type: 'text/plain'}), { ext:'pdb'}).then(function (o) {
        o.addRepresentation(m.representation, { color: scheme })
        o.autoView()
    });
    var scheme = NGL.ColormakerRegistry.addSelectionScheme(m.color_list, "new scheme");


    stage.setSpin(m.spin)
    document.addEventListener('spin', function(){
       stage.setSpin(m.spin);
    });

    window.addEventListener( "resize", function(){
        stage.handleResize();
    }, false );

    document.addEventListener('representation', function(){
        stage.compList[0].removeAllRepresentations();
        stage.compList[0].addRepresentation(m.representation, { color: scheme })
    });

    document.addEventListener('rcsb_id', function(){
        stage.removeAllComponents("");
        stage.loadFile(m.rcsb_id).then(function (o) {
            o.addRepresentation(m.representation, { color: scheme })
            o.autoView()
        });
    });

    document.addEventListener('color_list', function(){
        console.log(m.color_list)
        var list: Array<Array<String>> = m.color_list
        try{
              scheme = NGL.ColormakerRegistry.addSelectionScheme(list, "new scheme");
              stage.compList[0].reprList[0].setParameters( { color: scheme } );
        }
        catch(err) {
            console.log("badly defined color")
        }
    });

    document.addEventListener('pdb_string', function(){
        stage.removeAllComponents("");
        stage.loadFile( new Blob([m.pdb_string], {type: 'text/plain'}), { ext:'pdb'}).then(function (o) {
            o.addRepresentation(m.representation, { color: scheme })
            o.autoView()
        });

    });
   }

  get child_models(): LayoutDOM[] {
    return []
  }

  _update_layout(): void {
    this.layout = new LayoutItem()
    this.layout.set_sizing(this.box_sizing())
  }
}

export namespace ngl {
  export type Attrs = p.AttrsOf<Props>

  export type Props = LayoutDOM.Props & {
    spin: p.Property<boolean>
    representation: p.Property<string>
    rcsb_id: p.Property<string>
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

  static __name__ = "ngl"

  static init_ngl() {
    // This is usually boilerplate. In some cases there may not be a view.
    this.prototype.default_view = NGLView

    this.define<ngl.Props>(({String, Boolean, Any}) => ({
      spin:             [Boolean, false],
      representation:   [String],
      rcsb_id:          [String],
      color_list:       [Any],
      pdb_string:       [String]
    })
    )
  }
}
