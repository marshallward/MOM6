==================
MOM6 documentation
==================

MOM6 is the sixth generation of the Modular Ocean Model. This documentation
connects the model's scientific description and user guides with a browsable
reference to its Fortran implementation.

Start here
----------

* `About the documentation <about.html>`_
* `Working with MOM6 <working_with_MOM6.html>`_
* `Organization of the code <code_organization.html>`_: source directories,
  drivers, and build-time configuration.
* `Testing MOM6 <testing.html>`_
* `API reference <api-index.html>`_: modules, procedures, and call graphs.

Science and numerics
--------------------

* `Equations <equations.html>`_
* `Spatial discretization <discrete_space.html>`_
* `Time discretization <discrete_time.html>`_
* `Tracers <tracers.html>`_
* `Grids <grids.html>`_
* `Parameterizations <parameterizations.html>`_
* `Other physics <other_physics.html>`_
* `Forcing <forcing.html>`_
* `Parallel implementation <parallel.html>`_
* `Bibliography <zzbibliography.html>`_

About this site
---------------

This is an experimental documentation site. The handwritten manual chapters
and bibliography are available here. Detailed science pages generated from
``.dox`` files still open on Read the Docs; those links retain their original
page and section destinations during the migration.

The local API reference is generated from the source using Flint. Its
module-grouped caller and callee diagrams have expandable text views and
clickable procedure links. They show possible calls, not execution order or
the choices made by a particular runtime configuration.
