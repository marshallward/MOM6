Experimental Flint API deployment
=================================

The ``Flint API documentation`` workflow builds the reference independently of
the numerical and Doxygen checks. It scans the Fortran files under ``src`` using
Flint's default preprocessor configuration. It does not combine the alternative
drivers and infrastructure implementations in ``config_src``.

Build and deploy stages
-----------------------

Changes under ``src``, ``pkg``, or ``docs`` (or to ``.gitmodules``) trigger a build
on pushes and pull requests. Checkout includes submodules, which supply the
targets of source symlinks such as the TEOS10 files.
The workflow can also be started manually from the Actions page. Each build:

* Installs the Flint revision and Docutils version in ``requirements-api.txt``.
* Parses one file at a time, retaining documentation snapshots for linking.
* Generates one standard reST document per module and an index.
* Renders HTML with Docutils and fails on renderer warnings or broken page/anchor
  links. Unresolved symbols remain plain text and are recorded in diagnostics.
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

    git submodule update --init --recursive
    python3 -m venv .venv-api
    .venv-api/bin/python -m pip install -r docs/requirements-api.txt
    .venv-api/bin/python -m unittest discover -s docs -p 'test_build_api_docs.py'
    .venv-api/bin/python docs/build_api_docs.py --output docs/_build/flint-api

Open ``docs/_build/flint-api/html/index.html``. The output directory must be new
so stale pages cannot hide broken links; use another ``--output`` path for a
subsequent build, or remove the previous generated build directory yourself.

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

To preview theme changes without parsing Fortran again::

    .venv-api/bin/python docs/render_api_docs.py \
      docs/_build/flint-api/rst /tmp/mom6-api-preview

This changes presentation only. The generated reST is untouched and remains
usable with other renderers. The full build validates the theme's navigation
links as well as the API-reference links.
