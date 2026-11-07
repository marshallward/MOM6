Legacy documentation destination audit
======================================

Audit date: 2026-11-04.

Scope and result
----------------

The exact-pin release build contained ten distinct scientific/API destinations
on ``mom6.readthedocs.io``. Nine were inventory fallbacks from
``parameterizations_lateral.rst``; one was a hardcoded link in
``src/framework/Testing.rst``. All ten destinations have local content.

The migration initially mapped the nine audited Doxygen-qualified section
labels to their local targets while both builds coexisted. The chapters now
reference the source section labels directly, and that alias map has been
removed. The Testing chapter likewise uses a direct ``:doc:`` reference to
``mom_unit_scaling`` rather than a hardcoded published URL or generated path.

After these changes the full selected-source build has no scientific/API links
to Read the Docs. The redundant sidebar link to external "User documentation"
has also been removed now that the guide is local. The documentation home link
in ``about.rst`` still identifies Read the Docs as a hosting location.
The migration inventory has since been removed. Missing local destinations now
fail explicitly, including in partial builds with omitted API modules.

Content evidence
----------------

For each of the nine science destinations, the complete trailing namespace
comment block was extracted from the current Fortran source, converted using
the pinned Flint Doxygen frontend, and compared with generated reST. Every
converted block occurs intact in its owning module document: 976 comment lines
in total, containing 31 section/subsection commands. This checks the scientific
prose, equations, citations, lists, and subordinate sections rather than merely
the existence of an anchor. Docutils rendering and link validation also pass.

.. list-table:: Audited destinations
   :header-rows: 1

   * - Local page
     - Local source target
     - Namespace comment lines
   * - ``mom_hor_visc.html``
     - ``section_horizontal_viscosity``
     - 300
   * - ``mom_meke.html``
     - ``section_MEKE``
     - 229
   * - ``mom_mixed_layer_restrat.html``
     - ``section_mle``
     - 166
   * - ``mom_interface_filter.html``
     - ``section_interface_filter``
     - 34
   * - ``mom_lateral_mixing_coeffs.html``
     - ``section_Resolution_Function``
     - 79
   * - ``mom_tidal_forcing.html``
     - ``section_tides``
     - 32
   * - ``mom_load_love_numbers.html``
     - ``section_Love_numbers``
     - 34
   * - ``mom_self_attr_load.html``
     - ``section_SAL``
     - 40
   * - ``mom_spherical_harmonics.html``
     - ``section_spherical_harmonics``
     - 62
   * - ``mom_unit_scaling.html``
     - Module page (from Testing's See also link)
     - Not a namespace-section destination

A supplementary scan of non-symlink ``.F90`` files in the selected source roots
found 88 explicit section/subsection/anchor command occurrences. Every named
target appears in generated reST. This is a structural coverage check, not a
claim that every API declaration or unmarked source comment is documented.

Verification at the original audit
----------------------------------

* Nine builder/renderer tests pass, including local legacy-label resolution,
  owning-page validation, rejection of missing destinations, and root-relative API chapter
  links.
* Fresh exact-pin build: 278 source files, 278 module pages, 340 HTML pages,
  125,233 checked links, and zero unresolved call references.
* Graphs were disabled for this content audit; their layout is independent of
  reference resolution and scientific-comment conversion.
* The initial HTML scan left only two documentation home URLs on Read the Docs;
  the redundant sidebar link was subsequently removed as described above.

The remaining API completeness work (declaration details, undocumented
arguments, source links, and comments outside extracted documentation) is
separate from this legacy-destination audit. This audit compares current source
content with the new output; it does not claim byte-for-byte parity with an
older published documentation revision.

No Fortran code or Fortran comments were modified. Numerical answers are
unchanged.
