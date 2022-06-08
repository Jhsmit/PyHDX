"""
#TODO this file doesnt work from console, check paths
takes requirements from setup.cfg and makes req-<selection>.txt file
which can be used with `conda install --file req-<selection>.txt

selections are:
all
base
docs
pdf
web

"""
#%%

from __future__ import annotations


from configparser import ConfigParser
from functools import reduce
from operator import add
from pathlib import Path
from typing import Optional



#%%
# Pycharm scientific mode compat
if "__file__" not in locals():
    __file__ = Path().resolve() / "dev" / "deps" / "_requirements.py"
cwd = Path(__file__).parent
cwd

#%%

conversions = {"torch": "pytorch"}


def convert(req_list: list, reverse: bool = False):
    """
    Apply conversion dictionary. Default is pip -> conda (pytorch -> torch)
    """
    if reverse:
        c_dict = {v: k for k, v in conversions.items()}
    else:
        c_dict = conversions

    return [c_dict.get(item, item) for item in req_list]


#%%

EXTRAS = ["web", "pdf", "docs", "dev"]


def read_setup_cfg():
    cp = ConfigParser()
    setup_file = cwd.parent.parent / "setup.cfg"
    cp.read_string(setup_file.read_text())

    raw_deps = {}

    base = [p for p in cp.get("options", "install_requires").split("\n") if p]
    raw_deps["base"] = base

    for extra in EXTRAS:
        out = [p for p in cp.get("options.extras_require", extra).split("\n") if p]
        raw_deps[extra] = out

    return raw_deps


#%%


def remove_version_spec(dep_dict: dict) -> dict:
    """Removes the version specifier from dependency spec"""
    deps_out = {}
    for k, v in dep_dict.items():
        _deps = []
        for spec in v:
            # remove version spec part
            s = spec.replace("<", "").replace(">", "").split("=")[0]
            _deps.append(s.strip())
        deps_out[k] = _deps

    return deps_out


#%%


def conda_to_pip(conda_yml: Path, pip_txt: Path, extras: Optional[list[str]] = None):
    """Convert a conda .yml environement file to a pip requirements.txt file"""
    extras = extras or EXTRAS
    selection = ["base"] + extras

    raw_deps = read_setup_cfg()
    selected = {k: v for k, v in raw_deps.items() if k in selection}
    deps = remove_version_spec(selected)
    deps_flat = [p for d in deps.values() for p in d]

    import yaml

    linux_dict = yaml.safe_load((cwd / conda_yml).read_text())
    linux_split = [p.split("=") for p in linux_dict["dependencies"]]
    linux_deps, linux_versions, linux_id = zip(*linux_split)
    linux_pip = convert(linux_deps, reverse=True)

    overlap = set(linux_pip) & set(deps_flat)
    assert overlap == set(deps_flat)

    output = []
    for dep in deps_flat:
        idx = linux_pip.index(dep)
        version = linux_versions[idx]
        out = f"{dep}=={version}"
        output.append(out)

    (cwd / pip_txt).write_text("\n".join(output))


def make_requirements_files(extras: Optional[list[str]] = None):
    extras = extras or EXTRAS
    selection = ["base"] + extras
    raw_deps = read_setup_cfg()
    selected = {k: convert(v) for k, v in raw_deps.items() if k in selection}

    for s, dep_list in selected.items():
        (cwd / Path(f"req-{s}.txt")).write_text("\n".join(sorted(dep_list)))

    all_deps = set(reduce(add, selected.values()))
    (cwd / Path(f"req-all.txt")).write_text("\n".join(sorted(all_deps)))


#%%


if __name__ == "__main__":
    make_requirements_files()
    #
    platforms = ["linux", "windows"]
    for os in platforms:
        conda_file = Path(f"pinned/py38_{os}_conda.yml")
        pip_file = Path(f"pinned/py38_{os}_pip.txt")


        conda_to_pip(conda_file, pip_file)
