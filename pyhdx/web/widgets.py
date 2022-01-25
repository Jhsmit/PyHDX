import panel as pn
import param
from dask.distributed import as_completed
from panel.pane import HTML, Markdown
from panel.reactive import ReactiveHTML
from panel.util import as_unicode
from panel.widgets.input import _BkTextInput, StaticText


class NumericInput(pn.widgets.input.Widget):
    """
    NumericInput allows input of floats with bounds
    """

    type = param.ClassSelector(default=None, class_=(type, tuple), is_instance=True)

    value = param.Number(default=None)

    start = param.Number(default=None, allow_None=True)

    end = param.Number(default=None, allow_None=True)

    _rename = {
        "name": "title",
        "type": None,
        "serializer": None,
        "start": None,
        "end": None,
    }

    _source_transforms = {"value": """JSON.parse(value.replace(/'/g, '"'))"""}

    _target_transforms = {
        "value": """JSON.stringify(value).replace(/,/g, ", ").replace(/:/g, ": ")"""
    }

    _widget_type = _BkTextInput

    def __init__(self, **params):

        super(NumericInput, self).__init__(**params)
        self._state = ""
        self._validate(None)
        self._callbacks.append(self.param.watch(self._validate, "value"))

    def _validate(self, event):
        if self.type is None:
            return
        new = self.value
        if not isinstance(new, self.type) and new is not None:
            if event:
                self.value = event.old
            types = (
                repr(self.type) if isinstance(self.type, tuple) else self.type.__name__
            )
            raise ValueError(
                "LiteralInput expected %s type but value %s "
                "is of type %s." % (types, new, type(new).__name__)
            )

    def _bound_value(self, value):
        if self.start is not None:
            value = max(value, self.start)
        if self.end is not None:
            value = min(value, self.end)
        return value

    def _process_property_change(self, msg):
        if "value" in msg and msg["value"] is not None:
            try:
                value = float(msg["value"])
                msg["value"] = self._bound_value(value)
                if msg["value"] != value:
                    self.param.trigger("value")
            except ValueError:
                msg.pop("value")
        if "placeholder" in msg and msg["placeholder"] is not None:
            try:
                msg["placeholder"] = self._format_value(float(msg["placeholder"]))
            except ValueError:
                msg.pop("placeholder")
        return msg

    def _process_param_change(self, msg):
        msg = super(NumericInput, self)._process_param_change(msg)

        if "start" in msg:
            start = msg.pop("start")
            self.param.value.bounds[0] = start
        if "end" in msg:
            end = msg.pop("end")
            self.param.value.bounds[1] = end

        if "value" in msg:
            value = "" if msg["value"] is None else msg["value"]
            value = as_unicode(value)
            msg["value"] = value
        msg["title"] = self.name
        return msg


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


class NGL(ReactiveHTML):
    """
    The [NGL Viewer](https://github.com/nglviewer/ngl) can be used
    to show and analyse pdb molecule structures

    See also panel-chemistry for bokeh implementation:
    https://github.com/MarcSkovMadsen/panel-chemistry



    """

    object = param.String()

    extension = param.Selector(
        default="pdb",
        objects=EXTENSIONS,
    )

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
    completed = param.Integer(default=0)
    num_tasks = param.Integer(default=10, bounds=(1, None))

    async def run(self, futures):
        async for task in as_completed(futures):
            with pn.io.unlocked():
                self.completed += 1

    @property
    def value(self):
        value = int(100 * (self.completed / self.num_tasks))
        return max(
            0, min(value, 100)
        )  # todo check why this is somethings out of bounds

    def reset(self):
        self.completed = 0

    def increment(self):
        self.completed += 1

    @param.depends("completed", "num_tasks")
    def view(self):
        if self.value:
            return pn.widgets.Progress(
                active=True,
                value=self.value,
                align="center",
                sizing_mode="stretch_width",
            )
        else:
            return None


class CallbackProgress(pn.widgets.Progress):

    # this does not work as the worker is a different thread/process
    # see: https://distributed.readthedocs.io/en/latest/scheduling-state.html
    # and: https://stackoverflow.com/questions/44014988/how-to-get-information-about-a-particular-dask-task

    # also, when using ThreadPoolExecutor, the whole story becomes again different

    def callback(self, epoch, model, optimizer_obj):
        self.value = epoch
