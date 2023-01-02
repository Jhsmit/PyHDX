from functools import partial
from pathlib import Path
from urllib.request import urlretrieve

import yaml
from diagrams import Cluster, Diagram
from diagrams.aws.compute import Outposts, EC2, AutoScaling
from diagrams.aws.cost import CostExplorer
from diagrams.aws.database import Timestream, Redshift
from diagrams.aws.management import CommandLineInterface
from diagrams.custom import Custom

root_dir = Path(__file__).parent.parent.parent.parent

source_nodes = {"pyhdx": Timestream, "pdb": Redshift}

transform_nodes = {
    "table_source": Outposts,
    "cross_section": AutoScaling,
}
TRANSFORM_DEFAULT = EC2

opts_nodes = {}

OPTS_DEFAULT = CostExplorer

view_nodes = {
    "ngl_colors": partial(Custom, icon_path="molstar.png"),
    "logging": CommandLineInterface,
}

hv_logo_url = "https://holoviews.org/_static/logo.png"
hv_icon = "holoviews.png"
urlretrieve(hv_logo_url, hv_icon)
VIEW_DEFAULT = partial(Custom, icon_path=hv_icon)

sources = {}
transforms = {}
opts = {}
views = {}

cluster_attr = {
    "fontsize": "20",
    "bgcolor": "#e3e5e6",
}

diagram_attr = {
    # 'nodesep': '0.5',
    "layout": "dot",
    # 'pack': 'true',
    # 'clusterrank': 'local',
    # 'packMode': 'clust'
}


def make_diagram(name, yaml_dict):
    add_opts = False  # adding opts make the resulting scheme a bit messy
    opt_connections = []

    with Diagram(name, show=False, outformat="png", graph_attr=diagram_attr):
        for src_name, spec in yaml_dict.get("sources", None).items():
            node = source_nodes[spec["type"]]
            src = node(src_name)
            sources[src_name] = src

        for trs_name, spec in yaml_dict.get("transforms", None).items():
            node = transform_nodes.get(spec["type"], TRANSFORM_DEFAULT)
            trs = node(trs_name)
            transforms[trs_name] = trs
            if "source" in spec:
                source = sources.get(spec["source"]) or transforms.get(spec["source"])
                source >> trs
            if "sources" in spec:  # sources dict is label: src_name
                for src_id, src_name in spec["sources"].items():
                    source = sources.get(src_name) or transforms.get(src_name)
                    source >> trs

        # repeated code!
        for view_name, spec in yaml_dict.get("views", {}).items():
            node = view_nodes.get(spec["type"], VIEW_DEFAULT)
            view = node(view_name)
            views[view_name] = view
            if "source" in spec:
                source = sources.get(spec["source"]) or transforms.get(spec["source"])
                source >> view
            elif "sources" in spec:
                for src_id, src_name in spec["sources"].items():
                    source = sources.get(src_name) or transforms.get(src_name)
                    source >> view
            elif "views" in spec:
                for component_view_name in spec["views"]:
                    component_view = views.get(component_view_name)
                    component_view >> view
            if add_opts and "opts" in spec:
                for opt_name in spec["opts"]:
                    if isinstance(opt_name, dict):
                        pass
                    else:
                        try:
                            opt = opts[opt_name]
                            view << opt
                        except KeyError:
                            opt_connections.append((opt_name, view))

        for module_name, module_spec in yaml_dict["modules"].items():
            with Cluster(module_name, graph_attr=cluster_attr):
                for src_name, spec in module_spec.get("sources", {}).items():
                    node = source_nodes[spec["type"]]
                    src = node(src_name)
                    sources[src_name] = src

                for trs_name, spec in module_spec.get("transforms", {}).items():
                    node = transform_nodes.get(spec["type"], TRANSFORM_DEFAULT)
                    trs = node(trs_name)
                    transforms[trs_name] = trs
                    if "source" in spec:
                        source = sources.get(spec["source"]) or transforms.get(spec["source"])
                        source >> trs
                    if "sources" in spec:  # sources dict is label: src_name
                        for src_id, src_name in spec["sources"].items():
                            source = sources.get(src_name) or transforms.get(src_name)
                            source >> trs

                if add_opts:
                    for opt_name, spec in module_spec.get("opts", {}).items():
                        node = opts_nodes.get(spec["type"], OPTS_DEFAULT)
                        opt = node(opt_name)
                        opts[opt_name] = opt

                for view_name, spec in module_spec.get("views", {}).items():
                    node = view_nodes.get(spec["type"], VIEW_DEFAULT)
                    view = node(view_name)
                    views[view_name] = view
                    if "source" in spec:
                        source = sources.get(spec["source"]) or transforms.get(spec["source"])
                        source >> view
                    elif "sources" in spec:
                        for src_id, src_name in spec["sources"].items():
                            source = sources.get(src_name) or transforms.get(src_name)
                            source >> view
                    elif "views" in spec:
                        for component_view_name in spec["views"]:
                            component_view = views.get(component_view_name)
                            component_view >> view
                    if add_opts and "opts" in spec:
                        for opt_name in spec["opts"]:
                            if isinstance(opt_name, dict):
                                pass
                            else:
                                try:
                                    opt = opts[opt_name]
                                    view << opt
                                except KeyError:
                                    opt_connections.append((opt_name, view))

        if add_opts:
            with Cluster("Opts"):
                for opt_name, spec in yaml_dict["opts"].items():
                    node = opts_nodes.get(spec["type"], OPTS_DEFAULT)
                    opt = node(opt_name)
                    opts[opt_name] = opt

                for opt_name, view in opt_connections:
                    view << opts[opt_name]


if __name__ == "__main__":
    app_names = ["PyHDX_main_application", "PyHDX_rfu"]
    app_files = ["pyhdx_app.yaml", "rfu_app.yaml"]

    d = {}
    for name, file in zip(app_names, app_files):
        yaml_dir = root_dir / "pyhdx" / "web" / "apps" / file
        yaml_str = yaml_dir.read_text(encoding="utf-8")
        d[name] = yaml.safe_load(yaml_str)

    for name, yaml_dict in d.items():
        make_diagram(name, yaml_dict)
