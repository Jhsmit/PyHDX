import panel as pn
import param
from asyncio import as_completed
from panel.pane import HTML, Markdown
from panel.reactive import ReactiveHTML
from panel.widgets.input import StaticText


class ColoredStaticText(StaticText):
    _format = '<b>{title}</b>: <span class ="panel-colored-statictext">{value}</span>'


class HTMLTitle(HTML):
    title = param.String(doc="""Title""")
    priority = 0
    _rename = dict(HTML._rename, title=None)

    def __init__(self, **params):
        super().__init__(**params)
        self._update_title()

    @param.depends("title", watch=True)
    def _update_title(self):
        self.object = f"""<a class="title" href="" >{self.title}</a>"""


REPRESENTATIONS = [
    # "base",
    # "distance",
    "axes",
    "backbone",
    "ball+stick",
    "cartoon",
    "helixorient",
    "hyperball",
    "label",
    "licorice",
    "line",
    "point",
    "ribbon",
    "rocket",
    "rope",
    "spacefill",
    "surface",
    "trace",
    "unitcell",
    # "validation",
]
COLOR_SCHEMES = [
    "atomindex",
    "bfactor",
    "chainid",
    "chainindex",
    "chainname",
    "custom",
    "densityfit",
    "electrostatic",
    "element",
    "entityindex",
    "entitytype",
    "geoquality",
    "hydrophobicity",
    "modelindex",
    "moleculetype",
    "occupancy",
    "random",
    "residueindex",
    "resname",
    "sstruc",
    "uniform",
    "value",
    "volume",
]
EXTENSIONS = [
    "",
    "pdb",
    "cif",
    "csv",
    "ent",
    "gro",
    "json",
    "mcif",
    "mmcif",
    "mmtf",
    "mol2",
    "msgpack",
    "netcdf",
    "parm7",
    "pqr",
    "prmtop",
    "psf",
    "sd",
    "sdf",
    "top",
    "txt",
    "xml",
]

# todo remove
class NGL(ReactiveHTML):
    """
    The [NGL Viewer](https://github.com/nglviewer/ngl) can be used
    to show and analyse pdb molecule structures

    See also panel-chemistry for bokeh implementation:
    https://github.com/MarcSkovMadsen/panel-chemistry



    """

    object = param.String()

    extension = param.Selector(default="pdb", objects=EXTENSIONS,)

    background_color = param.Color(
        default="#F7F7F7", doc="Color to use for the background"
    )

    representation = param.Selector(
        default="ball+stick",
        objects=REPRESENTATIONS,
        doc="""
         A display representation. Default is 'ball+stick'. See
         http://nglviewer.org/ngl/api/manual/coloring.html#representations
         """,
    )

    color_scheme = param.Selector(default="chainid", objects=COLOR_SCHEMES)

    custom_color_scheme = param.List(
        default=[["#258fdb", "*"]],
        doc="""
        A custom color scheme. See
        http://nglviewer.org/ngl/api/manual/coloring.html#custom-coloring.""",
    )

    effect = param.Selector(
        default=None, objects=[None, "spin", "rock"], allow_None=True
    )

    _template = """
    <div id="ngl_stage" style="width:100%; height:100%;"></div>
    """
    _scripts = {
        "render": """
            var stage = new NGL.Stage(ngl_stage)        
            state._stage = stage
            state._stage.setParameters({ backgroundColor: data.background_color})
            stage.handleResize();
            self.updateStage()
        """,
        "object": """
            self.updateStage()
            """,
        "color_scheme": """
            self.setParameters()
            """,
        "custom_color_scheme": """
            self.setParameters()
        """,
        "background_color": """
        state._stage.setParameters({ backgroundColor: data.background_color})
        """,
        "setParameters": """
            if (state._stage.compList.length !== 0) {
                const parameters = self.getParameters();
                state._stage.compList[0].reprList[0].setParameters( parameters );
            }
            """,
        "getParameters": """
            if (data.color_scheme==="custom"){
                var scheme = NGL.ColormakerRegistry.addSelectionScheme( data.custom_color_scheme, "new scheme")
                var parameters = {color: scheme}
            }
            else {
                var parameters = {colorScheme: data.color_scheme}
            }
            parameters["sele"] = 'protein'

            return parameters
        """,
        "representation": """
            const parameters = self.getParameters();
            const component = state._stage.compList[0];
            component.removeAllRepresentations();
            component.addRepresentation(data.representation, parameters);
            """,
        "effect": """
            if (data.effect==="spin"){
                state._stage.setSpin(true);
            } else if (data.effect==="rock"){
                state._stage.setRock(true);
            } else {
                state._stage.setSpin(false);
                state._stage.setRock(false);
            }
            """,
        "updateStage": """
            parameters = self.getParameters();
            state._stage.removeAllComponents()
            state._stage.loadFile(new Blob([data.object], {type: 'text/plain'}), { ext: data.extension}).then(function (component) {
              component.addRepresentation(data.representation, parameters);
              component.autoView();
            });
            """,
        "after_layout": """
            state._stage.handleResize();
            """,
    }

    __javascript__ = [
        "https://unpkg.com/ngl@2.0.0-dev.38/dist/ngl.js",
    ]


class LoggingMarkdown(Markdown):
    def __init__(self, header, **params):
        super(LoggingMarkdown, self).__init__(**params)
        self.header = header
        self.contents = ""
        self.object = self.header + self.contents

    def write(self, line):
        self.contents = line + self.contents
        self.object = self.header + self.contents


class ASyncProgressBar(param.Parameterized):
    completed = param.Integer(default=0, doc="Number of completed jobs")

    num_tasks = param.Integer(default=10, doc="Total number of tasks", bounds=(1, None))

    active = param.Boolean(False, doc="Toggles the progress bar 'active' display mode")

    async def run(self, futures):
        self.active = True
        for task in as_completed(futures):
            await task
            self.active = False
            self.completed += 1

        self.reset()

    @property
    def value(self):
        value = int(100 * (self.completed / self.num_tasks))
        # todo check why this is sometimes out of bounds
        value = max(0, min(value, 100))

        if value == 0 and self.active:
            return None
        else:
            return value

    def reset(self):
        self.completed = 0

    def increment(self):
        self.completed += 1

    @param.depends("completed", "num_tasks", "active")
    def view(self):
        if self.value != 0:
            return pn.widgets.Progress(
                active=self.active,
                value=self.value,
                align="center",
                sizing_mode="stretch_width",
            )
        else:
            return pn.layout.Column()  # Or size 0 spacer?


class CallbackProgress(pn.widgets.Progress):

    # this does not work as the worker is a different thread/process
    # see: https://distributed.readthedocs.io/en/latest/scheduling-state.html
    # and: https://stackoverflow.com/questions/44014988/how-to-get-information-about-a-particular-dask-task

    # also, when using ThreadPoolExecutor, the whole story becomes again different

    def callback(self, epoch, model, optimizer_obj):
        self.value = epoch
