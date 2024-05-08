# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import os
from remix.framework import __version__


extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.extlinks",
    "sphinx.ext.ifconfig",
    "sphinx.ext.napoleon",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.jquery",
    "sphinx.ext.mathjax",
    "sphinx_multiversion",
]

bibtex_bibfiles = ["references.bib"]

myst_enable_extensions = [
    "dollarmath",
]
source_suffix = [".md", ".rst"]
master_doc = "index"
project = "REMix"
year = "2024"
author = "DLR"
copyright = "{0} {1}".format(author, year)
version = __version__

extlinks = {
    "issue": ("https://gitlab.com/dlr-ve/esy/remix/framework/-/issues/%s", "#%s"),
    "mr": ("https://gitlab.com/dlr-ve/esy/remix/framework/-/merge_requests/%s", "MR #%s"),
}

root_doc = "index"

pygments_style = "trac"

html_theme = "pydata_sphinx_theme"
html_favicon = "_static/images/DLR_Logo_REMix_short-name-only.svg"
html_sourcelink_suffix = ""

smartquotes = False
html_last_updated_fmt = "%b %d, %Y"
html_split_index = False

html_short_title = "%s-%s" % (project, version)
templates_path = ["_templates"]
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]


# Define the json_url for our version switcher.
if os.environ.get('CI_INSTANCE', '') == "dlr":
    json_url = "https://remix.pages.gitlab.dlr.de/framework/_static/versions.json"
else:
    json_url = "https://dlr-ve.gitlab.io/esy/remix/framework/_static/versions.json"
# 

# Set up the version switcher.  The versions.json is stored in the doc repo.
if os.environ.get('CIRCLE_JOB', False) and \
        os.environ.get('CIRCLE_BRANCH', '') != 'main':
    # For PR, name is set to its ref
    switcher_version = os.environ['CIRCLE_BRANCH']
elif ".dev" in version:
    switcher_version = "devdocs"
else:
    switcher_version = f"{version}"

html_theme_options = {
    "external_links": [
        {
            "url": "https://gitlab.com/dlr-ve/esy/remix/framework/-/releases",
            "name": "Changelog",
        },
    ],
    "gitlab_url": "https://gitlab.com/dlr-ve/esy/remix/framework/",
    "header_links_before_dropdown": 6,
    "logo": {
        "text": "",
        "image_light": "images/DLR_Logo_REMix_short-name-only.svg",
        "image_dark": "images/DLR_Logo_REMix_short-name-only.svg"
    },
    "use_edit_page_button": True,
    "show_toc_level": 4,
    "navbar_align": "left",
    "navbar_center": ["navbar-nav"],
    "announcement": "",
    "navigation_depth": 4,
    "show_nav_level": 4,
    "collapse_navigation": False,
    "navbar_start": ["navbar-logo"],
    "navbar_end": ["version-switcher", "theme-switcher", "navbar-icon-links"],
    "footer_items": ["copyright", "sphinx-version"],
    "switcher": {
        "version_match": switcher_version,
        "json_url": json_url,
    },
}

html_context = {
    "gitlab_user": "remix",
    "gitlab_repo": "framework",
    "gitlab_version": "dev",
    "doc_path": "docs",
}

napoleon_use_ivar = True
napoleon_use_rtype = False
napoleon_use_param = False


smv_tag_whitelist = r'^.*$'                          # Include all tags
smv_branch_whitelist = r'^(dev|docs\-test|pages.*)$' # Include dev branch or branches starting with pages
smv_remote_whitelist = r'^origin$'                   # Include any remote
smv_released_pattern = r'^\d+\.\d+.\d+$'             # All semantically versioned releases
