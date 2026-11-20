from pathlib import Path
import sys

DOCS_SOURCE_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = DOCS_SOURCE_DIR.parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT))

import packman

project = "PACKMAN"
copyright = "2026, PACKMAN contributors"
author = "PACKMAN contributors"
release = packman.__version__

root_doc = "index"
master_doc = root_doc

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_design",
]

autosummary_generate = False
autodoc_typehints = "description"
autoclass_content = "both"
napoleon_google_docstring = True
napoleon_numpy_docstring = True
todo_include_todos = False

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_static_path = ["../_static"]
html_css_files = ["custom.css"]
html_logo = "../_static/logo.svg"
html_favicon = "../_static/favicon.svg"
html_title = "PACKMAN documentation"
html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "navbar_align": "content",
    "navbar_start": ["navbar-logo"],
    "navbar_center": ["navbar-nav"],
    "navbar_end": ["theme-switcher"],
    "navigation_depth": 4,
    "show_nav_level": 2,
    "show_toc_level": 2,
}

copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True


def run_apidoc(_):
    from sphinx.ext import apidoc

    apidoc.main(
        [
            "-f",
            "-T",
            "-e",
            "-M",
            "-o",
            str(DOCS_SOURCE_DIR),
            str(REPOSITORY_ROOT / "packman"),
            str(REPOSITORY_ROOT / "packman" / "tests"),
        ]
    )


def setup(app):
    app.connect("builder-inited", run_apidoc)
