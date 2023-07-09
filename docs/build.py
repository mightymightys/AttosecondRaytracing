#!/usr/bin/env python3

import os
import shutil
import textwrap
from pathlib import Path

import pdoc.render_helpers

def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)

here = Path(__file__).parent

pdoc.render.configure(
    template_directory=here/ "pdoc" / "pdoc-template",
    math = True
)
# We can't configure Hugo, but we can configure pdoc.
pdoc.render_helpers.formatter.cssclass = "chroma pdoc-code"

modules = [
    #"DefaultOptions",
    "ModuleAnalysisAndPlots",
    "ModuleDetector",
    "ModuleGeometry",
    "ModuleMask",
    "ModuleMirror",
    "ModuleOpticalChain",
    "ModuleOpticalElement",
    "ModuleOpticalRay",
    "ModuleProcessing",
    "ModuleSource",
    "ModuleSupport",
    "ModuleDefects",
]
prefix = "../ART/"

# Generate API documentation using pdoc
api_output_dir = here / "pdoc" / "output"
if api_output_dir.exists():
    shutil.rmtree(api_output_dir)
api_output_dir.mkdir(parents=True)

pdoc.pdoc(*[prefix+m for m in modules], output_directory=api_output_dir)

# Generate markdown files from pdoc HTML files
api_content_dir = here / "src" / "content" / "API"
if api_content_dir.exists():
    shutil.rmtree(api_content_dir)
api_content_dir.mkdir(parents=True)
all_modules: dict[str, pdoc.doc.Module] = {}
for module_name in pdoc.extract.walk_specs([prefix+m for m in modules]):
    all_modules[module_name] = pdoc.doc.Module.from_name(module_name)

for module in modules:
    toc = pdoc.render.env.get_template("nav.html.jinja2").render(
        module=all_modules[pdoc.extract.walk_specs([prefix+module])[0]],
        all_modules=all_modules,
        root_module_name=pdoc.render_helpers.root_module_name(all_modules),
        )
    outfile = api_output_dir / "ART"/ f"{module.replace('.', '/')}_toc.html"
    outfile.parent.mkdir(parents=True, exist_ok=True)
    outfile.write_bytes(toc.encode())


    html_filename = api_output_dir / "ART" / f"{module}.html"
    markdown_filename = api_content_dir / f"{module}.md"
    markdown_content = textwrap.dedent(f"""\
        ---
        title: "{module}"
        url: "/api/{module}.html"
        math: true
        menu:
            addons:
                parent: 'API'
        ---

        {{{{< include-html "html/{module}.html" >}}}}
        """)
    markdown_filename.write_text(markdown_content)

API_index_content = """
---
title: "API"
linkTitle: "API"
weight: 4
---

{{% pageinfo %}}
The API may be updated without warning
{{% /pageinfo %}}

API stuff
"""
API_index_filename = api_content_dir / "_index.md"
API_index_filename.write_text(API_index_content)

site_output_dir = here / "src" / "content" / "API" / "html"
images_dir = here / "pdoc" / "images"
if site_output_dir.exists():
    shutil.rmtree(site_output_dir)

if (here / "src" / "static" / "api").exists():
    shutil.rmtree(here / "src" / "static" / "api")

htmls = api_output_dir / "ART"
os.mkdir(here / "src" / "static" / "api")
shutil.copytree(htmls, site_output_dir )
copytree(images_dir, here / "src" / "static" / "api" )
# Touch the addons-api.md file to trigger a rebuild of the addons API page
(here / "src" / "content" / "addons-api.md").touch()
