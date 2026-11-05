Flint documentation build and deployment
========================================

The ``Flint API documentation`` workflow builds the reference independently of
the numerical checks. It scans ``src`` together with these selected
configuration directories, using Flint's default preprocessor configuration:

* ``config_src/drivers/solo_driver``
* ``config_src/external`` (all subdirectories)
* ``config_src/infra/FMS2``
* ``config_src/memory/dynamic_symmetric`` (also supplies include headers)

Alternative drivers are not combined: for example, the solo and ice-solo drivers
both define ``atmos_ocean_fluxes_mod``. Duplicate module names fail the build
rather than overwriting pages or creating ambiguous links. FMS1 and asymmetric
memory are not included. External library implementations such as FMS itself
remain outside this source selection.

GSW and CVMix are external libraries, not MOM6 API pages. Files whose resolved
paths lie in ``pkg/GSW-Fortran`` or ``pkg/CVMix-src`` are excluded before parsing,
including TEOS10 symlinks under ``src``. MOM6's own EOS and CVMix wrappers remain
documented. ``reports/excluded-dependencies.txt`` records the excluded sources.

Build and deploy stages
-----------------------

The ``Linux verification`` workflow's ``check-style-and-docstrings`` job uses the
same pinned documentation tools instead of Doxygen/Sphinx. It runs the
builder/renderer tests and a complete content build with ``--graph-depth 0``,
including reference and HTML link validation. Graphviz is installed for the
small graph fixtures in the tests; the deployment workflow generates the full
set of diagrams. Check reports are retained as ``documentation-check-reports``
even when the content build fails. The job identifier is unchanged for existing
required-check configuration.

Changes under ``src``, ``config_src``, ``pkg``, or ``docs`` (or to ``.gitmodules``) trigger a build
on pushes and pull requests. The documentation checkout does not fetch
submodules. Symlinks into the excluded dependency trees are ignored even when
their targets are absent. Other broken source symlinks still fail the build.
The workflow can also be started manually from the Actions page. Each build:

* Installs the Flint, Docutils, and Pybtex dependencies pinned in
  ``requirements-api.txt``.
* Installs Graphviz and generates separate caller/callee diagrams, grouped by
  module, for connected procedures and interfaces.
* Parses one file at a time, retaining documentation snapshots for linking.
* Generates one standard reST document per module and an API index, and copies
  the maintained ``api-home.rst`` landing page and handwritten manual chapters.
* Renders API pages and manual chapters with Docutils, using the same site
  theme. Fails on renderer warnings or broken page/anchor
  links. External calls are grouped by dependency; unresolved symbols appear
  only in build diagnostics, not in published call lists.
* Uploads ``flint-api-reference`` containing ``rst``, ``html``, and ``reports``
  directories. PR builds can be downloaded for local preview.

The separate deployment job publishes successful builds from the repository's
default branch and the working branch ``check-assign-spacing``. It is restricted
to ``marshallward/MOM6``, which has no existing
Pages deployment at the time of setup. PRs and other feature branches only build
artifacts. Both deployment branches publish to the same site. The workflow does
not push generated files to a Git branch.

One-time GitHub setup
---------------------

In ``marshallward/MOM6``, open **Settings > Pages** and select **GitHub Actions**
as the build and deployment source. Ensure the ``github-pages`` environment
permits deployment from the default branch (currently ``dev/gfdl``) and
``check-assign-spacing``. If necessary, add the working branch under
**Settings > Environments > github-pages > Deployment branches and tags**.

Commit and push the workflow and build files on ``check-assign-spacing`` to test
deployment without merging into ``dev/gfdl``. A qualifying push will publish
the site at ``https://marshallward.github.io/MOM6/``. Manual dispatch is also
available once the workflow exists on the default branch. The deployment
job also reports the actual page URL. To adopt the workflow in another fork,
change the repository condition in ``.github/workflows/deploy-api.yml`` and
configure Pages there.

Read the Docs
-------------

``.readthedocs.yml`` uses custom build jobs to install the pinned dependencies
and Graphviz, run ``docs/build_api_docs.py``, and copy the resulting HTML into
``$READTHEDOCS_OUTPUT/html``. It does not invoke Sphinx or Doxygen. GitHub Pages
and Read the Docs therefore publish the same site structure.

The hosted build must be verified on Read the Docs after pushing this
configuration. Version activation and hosting settings remain managed there.
Older published versions keep their old builds. New versions use the standalone
site's page paths; redirecting old ``api/generated/...`` URLs is separate from
this build adapter.

Local build
-----------

Install Graphviz so that ``dot`` is on your PATH (for example,
``sudo apt-get install graphviz`` on Ubuntu). From the MOM6 repository root,
with Python 3.13::

    python3 -m venv .venv-api
    .venv-api/bin/python -m pip install -r docs/requirements-api.txt
    .venv-api/bin/python -m unittest discover -s docs -p 'test_build_api_docs.py'
    .venv-api/bin/python docs/build_api_docs.py --output docs/_build/flint-api

Open ``docs/_build/flint-api/html/index.html``. The output directory must be new
so stale pages cannot hide broken links; use another ``--output`` path for a
subsequent build, or remove the previous generated build directory yourself.

The home page introduces the documentation and links to the local science
and user-guide chapter pages. The local module index is
``api-index.html``; module page URLs and procedure anchors are unchanged.
Edit ``docs/api-home.rst`` to update the landing page. It uses ordinary reST.

``GUIDE_PAGES`` in ``render_guide_docs.py`` lists the 16 handwritten manual
chapters and their navigation titles, including the bibliography. They use
the maintained ``docs/*.rst`` sources, and appear under
"User and science guide" rather than in the API module list. Docutils renders
their HTML fragments for the existing MOM6 template. ``render_guide_docs.py``
resolves chapter navigation, named references, and citations across the site.

The 44 detailed science pages are maintained directly in ``src/**/*.rst``,
beside the implementation sources in their original ``.dox`` directories.
Their former ``src/**/*.dox`` sources have been retired. Edit the reST files;
these chapters require no markup conversion at build time.
``science-migration.json`` records the original paths, hashes, and conversion
counts. Each new path is the original source directory plus the recorded page
name and ``.rst`` extension. ``science_sources.py`` discovers these chapters
for the site builder. The migration preserved 493 inline formulae, 202 display equations,
35 figure placements, ten tables, and the scientific prose and citations.

Use standard ``:math:`` roles and ``.. math::`` directives for equations, with
explicit reST targets for equation references. Multiline equations use
``aligned``; keep their lines together without blank paragraphs. Reference links
use descriptive names (or the cited paper's equation numbers), rather than
unstable site-wide numbering. Figures use relative paths to ``docs/images/``
from their source directories; the builder adjusts these paths when staging.
The renderer downloads MathJax 3.2.2 (``tex-svg-full.js`` and its Apache-2.0
license) from version-pinned jsDelivr URLs and verifies their SHA-256 checksums.
Verified assets are cached in ``docs/_build/assets`` and copied into the site's
``assets`` directory. Local builds, GitHub Actions, and Read the Docs all use
this setup. The first build requires network access; subsequent builds can use
the verified cache offline. The generated site includes MathJax locally, so
equations render offline as well as online. URLs and checksums are maintained
in ``documentation_assets.py``; the third-party bundle is not stored in Git.

The new site requires neither Sphinx nor Doxygen. All reST is rendered with
Docutils. The existing ``toctree``, ``:ref:``, and ``:cite:`` source syntax is
supported by the standalone manual frontend. Chapter lists become ordinary links, with section links up
to the requested depth. Named references resolve against explicit Docutils
targets. Unknown references fail the build.

Fortran comments remain in their existing Doxygen format. The builder explicitly
selects Flint's ``doxygen`` comment frontend, which emits reST without running
Doxygen or Pandoc. Equations, tables, references, and literal diagrams are
interpreted during generation; source comments are not rewritten. Flint also
supports explicit ``plain`` and native ``rst`` modes, so a future source-format
migration is independent of replacing the build system. Unsupported Doxygen
commands fail generation with the owning source file in the error.

Generated API pages and manual chapters share the same citation and named-target
index. References to sections present in API comments resolve locally.
The legacy-destination audit in ``api-content-audit.rst`` records the content
checks and local mappings: the full selected-source build no longer uses the
old site for scientific/API destinations.

Pybtex reads the existing BibTeX files and formats the cited entries in the
``plain`` style. Citations link to stable, key-based bibliography anchors.
Missing citation keys and duplicate anchors fail the build. Pandoc and its
``citeproc`` filter are not required. The former Sphinx/Doxygen configuration,
extensions, build targets, and equation postprocessor have been retired.

The saved old-site inventory has been removed. Reference resolution needs no
network access and has no automatic external fallback. Unknown chapter
destinations and unresolved references fail the build. Partial source selections
must include the modules referenced by the maintained manual.

Repeat ``--source DIRECTORY`` to supply a custom set of source roots instead of
the defaults. For example, ``--source src`` restores a source-only reference.
The generated index records the selected directories.

Call graphs follow one hop by default; ``--graph-depth 2`` expands each
direction to two hops. Private and internal procedures are included so that
implementation paths remain visible. External and unresolved calls are omitted
from diagrams. These are possible-call relationships, not execution order.
DOT and SVG files are retained beside the reST pages. HTML embeds SVG inline
to preserve clickable node links, which are validated with the other site
links. ``reports/summary.json`` includes the number of generated diagrams.

Update the pinned Flint revision explicitly when adopting generator changes.
Generation and rendering remain separate functions in ``build_api_docs.py``;
the retained reST can be processed by another renderer. HTML destinations in
these particular files assume sibling ``.html`` pages.

This deployment changes documentation only; model numerical answers are
unchanged.

HTML presentation
-----------------

``render_api_docs.py`` assembles the rendered API and manual pages with the local
HTML template, stylesheet, and optional navigation script in ``api_theme``.
The theme provides a filterable module sidebar, breadcrumbs, in-page navigation,
procedure signature panels, responsive layouts, and print styling. Assets are
bundled with the site; there are no external fonts or JavaScript dependencies.
Body text, signatures, and graph labels use a 16px reading baseline; navigation
and secondary text use 14px. Graphviz labels are generated at 12pt (16px).
The ``--text-scale`` CSS variable scales both text and diagrams. Phone layouts
retain these reading sizes and let wide graphs scroll horizontally rather than
shrinking their labels. Graph scroll regions also accept keyboard focus.
Callee and caller diagrams each have their own initially collapsed disclosure
box. These work without JavaScript, and retain the scroll region when expanded.
Each graph box includes a collapsed ``Text view`` of the same caller-to-callee
edges, with qualified, clickable procedure names. The corresponding standalone
call lists are removed from HTML when a diagram is present; the original reST
retains them for other renderers. External-library calls remain separate.
Links and navigation work without JavaScript; the script adds module filtering
and initially collapses the module browser on small screens.

Procedure entries are grouped by standard reST containers. JavaScript turns
them into initially collapsed, keyboard-accessible disclosure panels. Following
a procedure link opens its panel (including enclosing panels for nested
procedures); printing opens all procedure panels temporarily. Without JavaScript,
the complete reference remains visible. Names retain definition spelling in
titles, navigation, signatures, and resolved links, while page filenames and
anchors keep their stable identities.

Top-level module sections are independently foldable through their headings.
They start open; individual procedure panels start closed. Following a link
opens every enclosing section and procedure panel. Printing temporarily opens
all content panels, then restores their previous state.

To preview theme changes without parsing Fortran again::

    .venv-api/bin/python docs/render_api_docs.py \
      docs/_build/flint-api/rst /tmp/mom6-api-preview

This changes presentation only. The generated reST is untouched and remains
usable with other renderers. The full build validates the theme's navigation
links as well as the API-reference links.

Dependency classification
-------------------------

``api-dependencies.toml`` specifies dependency implementation trees to exclude
from publication in its ``source_libraries`` table, relative to the repository
root. The file is optional; without it, GSW and CVMix implementation trees are
still excluded by default.

Explicit ``USE`` imports from unindexed modules identify external references
without requiring a module ownership inventory. Renames and re-exports through
documented modules retain the original module and symbol names. Here external
means outside the indexed source, not necessarily outside MOM6. This is import
provenance, not verification of a symbol's kind or a procedure's signature.
An optional ``libraries`` table can supply display labels, but MOM6 does not
maintain that mapping. The pinned Flint revision includes this
reference-classification behavior.
Wildcard imports without a symbol inventory and ambiguous references remain
unresolved; missing symbols are never classified by name prefix alone.

The public call lists contain resolved public project targets. External calls
appear in separate module groups (or library groups when labels are supplied),
collapsed when JavaScript is available.
They have no guessed documentation URLs. ``reports/external.txt`` records
their module and original symbol names; ``reports/unresolved.txt`` retains
the remaining unknowns for investigation. The build summary counts external
and unresolved references separately. Zero unresolved references remains the
goal; omitting unknown names from the site does not resolve them.
