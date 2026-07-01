import os
import sys

sys.path.insert(0, os.path.abspath("../../subprojects/doxyrest/sphinx"))

project = "rsx-rs"
copyright = "2024--present, rsx-rs developers"
author = "Rohit Goswami, Ruhila Goswami"
release = "0.2.5"

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
html_favicon = "_static/favicon.png"
html_title = "rsx-rs documentation"
html_baseurl = "https://rsx.rgoswami.me/"
# CSS via setup() so it loads after theme (readcon-core pattern)
# html_css_files left empty on purpose

# Edit-this-page + repo-stats (Shibuya / LODE readcon-core pattern).
# Narrative lives in docs/orgmode/; sidebar override points edit links at .org.
# Generated RST under docs/source/ stays gitignored (build artifact).
html_context = {
    "source_type": "github",
    "source_user": "HaoZeke",
    "source_repo": "rsx-rs",
    "source_version": "main",
    "source_docs_path": "/docs/orgmode/",
}

html_sidebars = {
    "**": [
        "sidebars/localtoc.html",
        "sidebars/repo-stats.html",
        "sidebars/edit-this-page.html",
    ],
}

html_theme_options = {
    "github_url": "https://github.com/HaoZeke/rsx-rs",
    "accent_color": "indigo",
    "dark_code": True,
    "globaltoc_expand_depth": 1,
    "light_logo": "_static/logo-light.svg",
    "dark_logo": "_static/logo-dark.svg",
    "og_image_url": "https://rsx.rgoswami.me/_static/og-image.png",
    "nav_links": [
        {"title": "Quickstart", "url": "tutorials/quickstart"},
        {
            "title": "Bindings",
            "children": [
                {
                    "title": "Language bindings",
                    "url": "reference/bindings",
                    "summary": "CLI, Python (pyrsx), R (rsxr), C matrix",
                },
                {
                    "title": "Rust API",
                    "url": "reference/rust-api",
                    "summary": "rsxcore crate reference",
                },
                {
                    "title": "C API",
                    "url": "reference/c-api",
                    "summary": "cbindgen / Doxyrest",
                },
                {
                    "title": "R integration",
                    "url": "howto/r-integration",
                    "summary": "In-process rsxr and CLI subprocess",
                },
            ],
        },
        {
            "title": "Reference",
            "children": [
                {"title": "Commands", "url": "reference/commands", "summary": "CLI flags"},
                {"title": "Glossary", "url": "reference/glossary", "summary": "Terms"},
                {"title": "FAQ", "url": "howto/faq", "summary": "Common questions"},
            ],
        },
        {
            "title": "Explain",
            "children": [
                {"title": "Architecture", "url": "explanation/architecture"},
                {"title": "Benchmarks", "url": "explanation/benchmarks"},
                {"title": "Citation", "url": "explanation/citation"},
            ],
        },
    ],
}

copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: |R> "
copybutton_prompt_is_regexp = True
copybutton_exclude = ".linenos, .gp, .go"
copybutton_line_continuation_character = "\\"
copybutton_here_doc_delimiter = "EOF"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

rust_crates = {
    "rsx_core": os.path.abspath("../../rsxcore/"),
}
rust_doc_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "crates")
rust_rustdoc_fmt = "rst"
rust_generate_mode = "always"

rustdoc_postprocess_toctree_target = "reference/rust-api.rst"
rustdoc_postprocess_toctree_rst = """
Rust API (``rsx_core``)
-----------------------

.. toctree::
   :maxdepth: 2

   ../crates/rsx_core/lib
"""


def setup(app):
    # priority > default so we load *after* doxyrest-pygments.css (which paints
    # many Pygments names #000 !important and wins on equal/earlier cascade).
    app.add_css_file("custom.css", priority=1000)
