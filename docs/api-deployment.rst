Experimental Flint API deployment
=================================

The ``Flint API documentation`` workflow builds the reference independently of
the numerical and Doxygen checks. It scans ``src`` together with these selected
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

Changes under ``src``, ``config_src``, ``pkg``, or ``docs`` (or to ``.gitmodules``) trigger a build
on pushes and pull requests. The documentation checkout does not fetch
submodules. Symlinks into the excluded dependency trees are ignored even when
their targets are absent. Other broken source symlinks still fail the build.
The workflow can also be started manually from the Actions page. Each build:

* Installs the Flint revision and Docutils version in ``requirements-api.txt``.
* Parses one file at a time, retaining documentation snapshots for linking.
* Generates one standard reST document per module and an index.
* Renders HTML with Docutils and fails on renderer warnings or broken page/anchor
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

Local build
-----------

From the MOM6 repository root, with Python 3.13::

    python3 -m venv .venv-api
    .venv-api/bin/python -m pip install -r docs/requirements-api.txt
    .venv-api/bin/python -m unittest discover -s docs -p 'test_build_api_docs.py'
    .venv-api/bin/python docs/build_api_docs.py --output docs/_build/flint-api

Open ``docs/_build/flint-api/html/index.html``. The output directory must be new
so stale pages cannot hide broken links; use another ``--output`` path for a
subsequent build, or remove the previous generated build directory yourself.

Repeat ``--source DIRECTORY`` to supply a custom set of source roots instead of
the defaults. For example, ``--source src`` restores a source-only reference.
The generated index records the selected directories.

Update the pinned Flint revision explicitly when adopting generator changes.
Generation and rendering remain separate functions in ``build_api_docs.py``;
the retained reST can be processed by another renderer. HTML destinations in
these particular files assume sibling ``.html`` pages.

This deployment changes documentation only; model numerical answers are
unchanged.

HTML presentation
-----------------

``render_api_docs.py`` renders the generated reST using Docutils and the local
HTML template, stylesheet, and optional navigation script in ``api_theme``.
The theme provides a filterable module sidebar, breadcrumbs, in-page navigation,
procedure signature panels, responsive layouts, and print styling. Assets are
bundled with the site; there are no external fonts or JavaScript dependencies.
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

``api-dependencies.toml`` records the ownership of known dependency modules.
The file is optional. Without it, the builder still excludes the GSW and CVMix
implementation trees by default, but imports from unindexed modules remain
in unresolved diagnostics rather than receiving guessed library labels.
Explicit imports from unindexed modules are identified as potentially external
in those diagnostics, with their original module and symbol names. The build
prints a warning and counts them as ``potential_external_references`` (a subset
of ``unresolved_references``). An unindexed module may also be an omitted project
source, so it is not automatically classified as a confirmed dependency.
With older Flint revisions that lack import-provenance metadata, this count is
``null`` and the builder issues a general warning when unresolved calls exist.
The ``libraries`` and ``source_libraries`` tables are independently optional.
Its ``source_libraries`` table also specifies dependency implementation trees
to exclude from publication, relative to the MOM6 repository root.
An explicit ``USE`` import from one of these modules identifies an external
reference, preserving renames and re-exports through documented modules. This
is import-based provenance, not verification of the dependency's interface.
Wildcard imports without a symbol inventory and ambiguous references remain
unresolved; missing symbols are never classified by name prefix alone.

The public call lists contain resolved public project targets. External calls
appear in separate library groups, collapsed when JavaScript is available.
They have no guessed documentation URLs. ``reports/external.txt`` records
their module and original symbol names; ``reports/unresolved.txt`` retains
the remaining unknowns for investigation. The build summary counts external
and unresolved references separately. Zero unresolved references remains the
goal; omitting unknown names from the site does not resolve them.
