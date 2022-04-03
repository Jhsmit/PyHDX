import collections

from distributed import Client

from pyhdx.local_cluster import default_client
from pyhdx.support import gen_subclasses
from pyhdx.web.controllers import *
from pyhdx.web.main_controllers import MainController
from pyhdx.web.opts import OptsBase
from pyhdx.web.sources import *
from pyhdx.web.tools import supported_tools
from pyhdx.web.transforms import *
from pyhdx.web.views import View
from pyhdx.web.cache import Cache

element_count = 0


class AppConstructor(param.Parameterized):

    sources = param.Dict(default={})

    transforms = param.Dict(default={})

    opts = param.Dict(default={})

    tools = param.Dict(default={})

    views = param.Dict(default={})

    loggers = param.Dict(default={})

    ctrl_class = param.ClassSelector(class_=MainController, instantiate=False)

    client = param.ClassSelector(default=None, class_=Client)

    cache = param.ClassSelector(default=Cache(), class_=Cache)

    def __init__(self, **params):
        super().__init__(**params)
        self.classes = self.find_classes()

    def parse(self, yaml_dict, **kwargs):
        self._parse_sections(yaml_dict)
        for name, dic in yaml_dict.get("modules", {}).items():
            self._parse_sections(dic)

        d = yaml_dict["controllers"]
        controllers = {name: self._resolve_class(name, "controller") for name in d}

        main_ctrl = yaml_dict["main_controller"]
        _type = main_ctrl.pop("type")
        main_ctrl_class = self._resolve_class(_type, "main")
        ctrl = main_ctrl_class(
            controllers.values(),
            sources=self.sources,
            transforms=self.transforms,
            opts=self.opts,
            views=self.views,
            loggers=self.loggers,
            **kwargs,
            **main_ctrl,
        )

        return ctrl

    @staticmethod
    def find_classes():
        base_classes = {
            "main": MainController,
            "transform": Transform,
            "source": Source,
            "view": View,
            "opt": OptsBase,
            "controller": ControlPanel,
        }
        classes = {}
        for key, cls in base_classes.items():
            base_cls = base_classes[key]
            all_classes = list(
                [cls for cls in gen_subclasses(base_cls) if getattr(cls, "_type", None)]
            )
            all_classes.append(base_cls)
            types = [cls._type for cls in all_classes]
            if len(types) != len(set(types)):
                duplicate_items = [
                    item
                    for item, count in collections.Counter(types).items()
                    if count > 1
                ]
                raise ValueError(
                    f"Multiple implementations of {key!r} found with the same type: {duplicate_items}"
                )
            class_dict = {cls._type: cls for cls in all_classes}
            classes[key] = class_dict

        classes["tool"] = supported_tools

        return classes

    def _parse_sections(self, yaml_dict):
        sections = ["sources", "transforms", "tools", "opts", "views"]
        for section in sections:
            element = section[:-1]
            element_dict = getattr(self, element + "s")

            d = yaml_dict.get(section, {})
            for name, spec in d.items():
                if name in element_dict:
                    raise ValueError(
                        f"The element {element!r} with name {name!r} already exists"
                    )
                # todo move to classmethod on object which checks spec/kwargs  (also prevents logger from needing a source)
                if "type" not in spec:
                    raise KeyError(
                        f"The field 'type' is not specified for {section[:-1]} {name!r}"
                    )
                # _type = spec.pop('type')
                if section in ["transforms", "views"] and "source" not in spec:
                    # raise KeyError(f"The field 'source' is not specified for {section[:-1]} {name!r}")
                    print(
                        f"The field 'source' is not specified for {section[:-1]} {name!r}"
                    )
                obj = self.create_element(name, element, **spec)
                element_dict[name] = obj

    def create_element(self, name: str, element: str, **spec):
        """

        :param name:
        :param element: eiter source, filter, opt, view, tool
        :param spec:
        :return:
        """
        global element_count

        _type = spec.pop("type")
        kwargs = self._resolve_kwargs(**spec)
        class_ = self._resolve_class(_type, element)
        if element == "transform":
            kwargs["_cache"] = self.cache
        obj = class_(name=name, **kwargs)
        element_count += 1

        return obj
        #

    def _resolve_class(self, _type, cls):
        return self.classes[cls][_type]

    def _resolve_kwargs(self, **kwargs):
        global element_count

        resolved = {}
        for k, v in kwargs.items():
            if k == "source":
                # temporary:
                if v is None:
                    resolved[k] = v
                else:
                    obj = self.sources.get(v) or self.transforms.get(
                        v
                    )  # can be none in case of logging
                    resolved[k] = obj
            elif k == "sources":
                # v should be a dict: src_type (view spec): src_name
                sources = {}
                for src_type, src in v.items():
                    obj = self.sources.get(src) or self.transforms.get(src)
                    sources[src_type] = obj
                # obj = {src_type: self.sources[src] for src_type, src in v.items()}
                resolved[k] = sources
            elif k == "opts":
                v = (
                    [v] if isinstance(v, (str, dict)) else v
                )  # allow single opt by str/dict (needs testing)
                opts = []
                for vi in v:
                    if isinstance(vi, dict):  # in situ opt declaration
                        if len(vi) != 1:
                            raise ValueError("Opts ")
                        name = next(iter(vi))  # get the first key
                        obj = self.create_element(
                            f"{name}_{element_count:05d}", "opt", **vi[name]
                        )  # should these in situ opts be added to global opts? probably not or they should have nested names
                        opts.append(obj)
                    else:
                        opts.append(self.opts[vi])

                resolved[k] = opts  # [self.opts[vi] for vi in v]
            elif k == "views":
                v = [v] if isinstance(v, str) else v  # allow single view by str
                resolved[k] = [self.views[vi] for vi in v]
            elif k == "tools":
                v = [v] if isinstance(v, str) else v  # allow single tool by str
                resolved[k] = [self.tools[vi] for vi in v]
            elif (
                k == "dependencies"
            ):  # dependencies are opts/transforms/controllers? (anything with .updated event)
                all_objects = []
                for type_, obj_list in v.items():
                    for obj in obj_list:
                        all_objects.append(getattr(self, type_)[obj])
                resolved[k] = all_objects
            elif k == "logger":
                resolved[k] = self.loggers[v]
            elif k == "tooltips":
                # workaround for pyyaml not reading tuples directly
                resolved[k] = [tuple(item) for item in v]

            else:
                resolved[k] = v

        return resolved
