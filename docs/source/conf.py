project = "rsx-rs"
copyright = "2024--present, rsx-rs developers"
author = "Rohit Goswami"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinxcontrib_rust",
    "sphinx_rustdoc_postprocess",
    "sphinx_click",
]

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "shibuya"
html_static_path = []

html_theme_options = {
    "github_url": "https://github.com/HaoZeke/rsx-rs",
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

rust_crates = {
    "rsx_core": "rsx-core/",
}
rust_doc_dir = "api/rust"
