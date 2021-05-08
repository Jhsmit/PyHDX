import panel as pn
import param
from awesome_panel_extensions.web_component import WebComponent

MWC_ICONS = [None, "accessibility", "code", "favorite"]

class NGLView(WebComponent):
    html = param.String("""
        <div id="viewport" style="width:100%; height:100%;"></div>
            <script>
            stage = new NGL.Stage("viewport");
            stage.loadFile("rcsb://1CRN").then(function(o){{ 
                o.addRepresentation("cartoon");
                o.autoView();
            }});            
            </script>
        """)
    attributes_to_watch = param.Dict({})
    events_to_watch = param.Dict({})


nglview = NGLView()
#mwc_button = MWCButton(name="Click Me!")

MWC_EXTENSIONS = """
            <script type='module' src='https://cdn.jsdelivr.net/gh/arose/ngl@v2.0.0-dev.37/dist/ngl.js'></script>
            """
extensions_pane = pn.pane.Markdown(MWC_EXTENSIONS,height=0, width=0,sizing_mode="fixed", margin= 0)

app = pn.Column(
    extensions_pane, nglview
)

app.servable()
