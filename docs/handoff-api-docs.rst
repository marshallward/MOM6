API documentation development handoff
=====================================

Updated 2026-11-03. This is a temporary development handoff; these development
commits are not intended to be preserved as the final contribution history.

Goal and current state
----------------------

Develop Flint as a general Fortran API documentation generator and use it to
publish a useful MOM6 reference. Flint generates portable reStructuredText;
MOM6 owns source selection, dependency metadata, HTML presentation, and deployment.
The HTML renderer uses Docutils and local CSS/JavaScript, without a frontend
framework or external runtime assets.

Repositories on the original computer:

* MOM6: ``/home/marshall/gfdl/mom6``, branch ``check-assign-spacing``, push remote
  ``marshall``. HEAD at handoff: ``451c52949`` (section folding and larger fonts).
* Flint: ``/home/marshall/python/flint``, branch ``main``, push remote ``marshall``.
  HEAD at handoff: ``49ca238c50069f9a65c2b9455422608bcd642729``.
* ``docs/requirements-api.txt`` pins that Flint commit and Docutils 0.22.2.
  Verify commits have been pushed before moving; local HEAD is not proof of that.

The only tracked uncommitted change before adding this handoff was in
``docs/api_theme/api.css``: restore ``--text-scale: 1`` and sidebar width 278px.
Flint had no tracked modifications. Both repositories have unrelated untracked
files: stage explicit paths, not ``git add .``.

Latest user decision
--------------------

The user liked independent section folding. They initially requested text about
50% larger, then decided that increase was premature and requested undoing it
everywhere. The CSS rollback restores the original 15px body text, 13px procedure
signatures, and 278px sidebar. Folding remains. The scale variable is retained
at 1; no mobile-only override was requested or added.

Typography should be revisited later, not redesigned now. The assistant initially
misinterpreted this as a request for a broader redesign and then over-explained
mobile overrides. Clarify intent before turning design discussion into edits.

Move and reproduce
------------------

Commit/push the intended work and clone both repositories on the new computer,
or otherwise transfer those changes explicitly. Recreate environments rather
than copying virtual environments. From the MOM6 checkout::

    python3 -m venv .venv-api
    .venv-api/bin/python -m pip install -r docs/requirements-api.txt
    .venv-api/bin/python -m unittest discover -s docs -p 'test_build_api_docs.py'
    .venv-api/bin/python docs/build_api_docs.py --output /tmp/mom6-api-preview

The output directory must not already exist. Use a new directory for each full
build. Open ``/tmp/mom6-api-preview/html/mom.html`` in a browser, or serve it::

    .venv-api/bin/python -m http.server 8000 --bind 127.0.0.1 --directory /tmp/mom6-api-preview/html

Then visit ``http://127.0.0.1:8000/mom.html``. Temporary previews from the old
computer are disposable and need not be transferred.

For active Flint development, substitute the new local checkout path::

    PYTHONPATH=/path/to/flint .venv-api/bin/python docs/build_api_docs.py --output /tmp/mom6-api-local

Always validate deployment using the exact installed requirements as well.
Testing only against local Flint previously hid a CI version mismatch.

For presentation-only changes, rerender existing reST without reparsing Fortran::

    .venv-api/bin/python docs/render_api_docs.py /tmp/mom6-api-preview/rst /tmp/mom6-api-restyled/html
    node --check docs/api_theme/api.js
    git diff --check

Important implementation files
------------------------------

* ``docs/build_api_docs.py``: selected sources, dependency exclusions, generation,
  diagnostics, and HTML link validation.
* ``docs/render_api_docs.py``: Docutils renderer, title-case labels, navigation.
* ``docs/api_theme/api.css`` and ``api.js``: typography, layout, folding,
  deep-link reveal, print-state restoration, and module filtering.
* ``docs/test_build_api_docs.py``: six builder/renderer tests.
* ``docs/api-dependencies.toml``: optional module-to-library ownership metadata.
* ``docs/api-deployment.rst``: detailed publishing and dependency policy.
* ``.github/workflows/deploy-api.yml``: GitHub Pages build/deploy.
* Flint ``flint/docgen/``: documentation model, resolution, destinations, reST.
* Flint ``flint/macros.py`` and ``flint/lexer.py``: preprocessing support.

Behavior to preserve
--------------------

* Display original case, preferring declaration spelling for resolved names.
  Keep stable lowercase filenames/anchors and compact ``func(a, b, c)`` signatures.
* Top-level module sections start open and fold independently. Procedure entries
  start closed. Native ``details``/``summary`` elements implement folding; no
  animation is needed. Without JavaScript, ordinary reST content stays visible.
* Direct links open enclosing sections/procedures, including repeated same-hash
  links. Printing temporarily expands content and restores the previous state.
* Hide private implementation targets from forward call/implementation lists;
  retain public generic-interface links. Explicit procedure interfaces, including
  ``bind(c)``, are documented even without bodies. Dummy procedures are arguments,
  not published call targets.
* Exclude GSW/CVMix implementations, including TEOS10 symlinks, but retain MOM6
  wrappers. Classify external ownership through actual USE associations and
  configured dependency modules, not symbol-prefix guesses.
* External calls have separate initially collapsed library groups. Unknown calls
  remain in build diagnostics rather than being assigned guessed targets.
* Optional missing dependency metadata retains default dependency exclusions.
  Older Flint without ``Reference.unindexed_module`` is supported and reports
  ``potential_external_references: null``.

Source selection and deployment
-------------------------------

Selected roots are ``src``, ``config_src/drivers/solo_driver``, all of
``config_src/external``, ``config_src/infra/FMS2``, and
``config_src/memory/dynamic_symmetric``. Do not combine drivers yet: solo and
ice-solo both define ``atmos_ocean_fluxes_mod``.

Documentation CI does not recursively fetch submodules. Known dependency
symlinks may be broken and must be excluded before parsing.

Published URL: ``https://www.marshallward.org/MOM6/index.html``.
GitHub Pages paths are case-sensitive: retain canonical ``/MOM6/``.

When Flint changes, publish Flint first, update the exact MOM6 requirements pin,
then test with a fresh installation of that pin. Previously omitted untracked
macro files caused a deployment failure, so verify required files are tracked.
The user has been handling commits/pushes, including interactive signing.

Validation history and known limitations
----------------------------------------

The last full reference build validated 278 selected source files, 279 HTML pages
including the index, 100,701 links, zero unresolved references, and 528 external
references (FMS 406, NetCDF 79, GSW 16, CVMix 27). It excluded 32 GSW source paths.

Before this handoff, six MOM6 builder tests, JavaScript syntax, and diff whitespace
checks passed. Browser checks passed for section/procedure folding, nested anchor
reveal, print-state restoration, and absence of horizontal overflow on the tested
desktop/mobile viewports. Those layout checks used the enlarged text. The final
font rollback has only had a diff whitespace check.

Flint's last broad run passed 268 runnable tests. Full discovery has an existing
unrelated failure in ``test/test_real_units.py`` importing missing
``report_real_units``; that file was excluded from the broad run.

Function-like preprocessing supports nested arguments, continuation lines,
chained expansion, recursion protection, conditionals, and undefinition while
preserving source round-trips. Variadic macros, stringification, token pasting,
and function-macro expansion in ``#if`` expressions are unsupported. Parentheses
and strings protect argument commas; square brackets do not.

Assignment-spacing performance work is deferred. A macro-free fast path restored
the measured run to about 17.3 seconds (from 18.1 seconds); the user deferred
investigation of the remaining roughly five-second overhead.

Next steps
----------

1. Transfer the intended commits and CSS rollback, and recreate the environment.
2. Run the exact-pin builder tests and generate a fresh preview.
3. Resume typography discussion only when requested; the current decision is
   original-size text with section/procedure folding.

Read repository agent instructions when resuming, including the MOM6 AI policy
reminder and code-style requirements. These documentation-only changes do not
alter numerical results; all answers are bitwise identical.
