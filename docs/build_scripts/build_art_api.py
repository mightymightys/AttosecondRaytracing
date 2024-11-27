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


def remove_module_name_div(html_path):
    with open(html_path, 'r', encoding='utf-8') as file:
        content = file.read()

    # Find and remove the modulename div
    start_marker = '<h1 class="modulename">'
    end_marker = '</h1>'
    start_index = content.find(start_marker)
    end_index = content.find(end_marker, start_index)
    
    if start_index != -1 and end_index != -1:
        content = content[:start_index] + content[end_index + len(end_marker):]

    with open(html_path, 'w', encoding='utf-8') as file:
        file.write(content)

here = Path(__file__).parent.parent

pdoc.render.configure(
    template_directory=here / "pdoc" / "pdoc-template",
    math=True,
    search=False,
    mermaid=True
)
# We can't configure Hugo, but we can configure pdoc.
pdoc.render_helpers.formatter.cssclass = "chroma pdoc-code"

modules = [
    "ARTmain",
    "ModuleAnalysisAndPlots",
    "ModulePlottingMethods",
    "ModuleAnalysis"
]
prefix = "../AttosecondRayTracing/src/ART/"

# Generate API documentation using pdoc
api_output_dir = here / "pdoc" / "output"
if api_output_dir.exists():
    shutil.rmtree(api_output_dir)
api_output_dir.mkdir(parents=True)

pdoc.pdoc(*[prefix+m for m in modules], output_directory=api_output_dir)

# Process each HTML file to remove the modulename div
for module in modules:
    html_path = api_output_dir / "ART" / f"{module}.html"
    remove_module_name_div(html_path)


# Generate markdown files from pdoc HTML files
api_content_dir = here / "src" / "content" / "ARTAPI"
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
        url: "/artapi/{module}.html"
        type: docs
        math: true
        menu:
            addons:
                parent: 'ART API'
        ---

        {{{{< include-html ".html/{module}.html" >}}}}
        """)
    markdown_filename.write_text(markdown_content)

API_index_content = """
---
title: "ART API"
linkTitle: "ART API"
weight: 4
type: docs
description: >
    User-facing package, contains mostly visualisation and analysis tools.
---

{{% pageinfo %}}
The API may be updated without warning
{{% /pageinfo %}}

API stuff TODO
"""
API_index_filename = api_content_dir / "_index.md"
API_index_filename.write_text(API_index_content)

site_output_dir = here / "src" / "content" / "ARTAPI" / ".html"
images_dir = here / "pdoc" / "images"
if site_output_dir.exists():
    shutil.rmtree(site_output_dir)

if (here / "src" / "static" / "api").exists():
    shutil.rmtree(here / "src" / "static" / "api")

htmls = api_output_dir / "ART"
os.mkdir(here / "src" / "static" / "api")
shutil.copytree(htmls, site_output_dir )
copytree(images_dir, here / "src" / "static" / "api" )
