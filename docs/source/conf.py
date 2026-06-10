import os
import sys

sys.path.insert(0, os.path.abspath("../../subprojects/doxyrest/sphinx"))

project = "rsx-rs"
copyright = "2024--present, rsx-rs developers"
author = "Rohit Goswami"
release = "0.2.3"

extensions = [
    "doxyrest",
    "cpplexer",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinxcontrib_rust",
    "sphinx_rustdoc_postprocess",
    "sphinx_click",
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx_tabs.tabs",
]

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "shibuya"
html_static_path = ["_static"]

html_theme_options = {
    "github_url": "https://github.com/HaoZeke/rsx-rs",
    "light_logo": "_static/logo-light.svg",
    "dark_logo": "_static/logo-dark.svg",
    "og_image_url": "https://rsx.rgoswami.me/_static/og-image.png",
}

html_favicon = "_static/favicon.png"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

import os

rust_crates = {
    "rsx_core": os.path.abspath("../../radsex-core/"),
}
rust_doc_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "crates")
rust_rustdoc_fmt = "rst"
rust_generate_mode = "always"

# -- sphinx-rustdoc-postprocess configuration --------------------------------
rustdoc_postprocess_toctree_target = "reference/rust-api.rst"
rustdoc_postprocess_toctree_rst = """
Rust API (``rsx_core``)
-----------------------

.. toctree::
   :maxdepth: 2

   ../crates/rsx_core/lib
"""

