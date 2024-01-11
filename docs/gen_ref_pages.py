"""Generate the code reference pages."""
from pathlib import Path

import mkdocs_gen_files

ROOT_DIR = "pyhdx"
nav = mkdocs_gen_files.Nav()

# open_func = open # for debugging
open_func = mkdocs_gen_files.open

# %%
skip = ["pyhdx_diagram"]
for path in sorted(Path(ROOT_DIR).rglob("*.py")):  #
    module_path = path.relative_to(ROOT_DIR).with_suffix("")  #
    doc_path = path.relative_to(ROOT_DIR).with_suffix(".md")  #

    full_doc_path = Path("reference", doc_path)  #

    parts = list(module_path.parts)
    if module_path.stem.startswith("_"):
        continue
    if module_path.stem in skip:
        continue

    nav[parts] = doc_path.as_posix()

    with open_func(full_doc_path, "w") as fd:  #
        identifier = ".".join(parts)  #
        print("::: " + identifier, file=fd)  #

    mkdocs_gen_files.set_edit_path(full_doc_path, path)  #


with open_func("reference/SUMMARY.md", "w") as nav_file:  #
    nav_file.writelines(nav.build_literate_nav())
