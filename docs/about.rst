About this documentation
========================

MOM6 documentation combines a Flint-generated API reference with the user and
science guide. The same standalone HTML site can be published on GitHub Pages
and Read the Docs.

Here is where to find particular documentation:

Download, compile and run
  Installation documentation is in the form of user-driven (editable) wiki attached to the MOM6-examples GitHub repository.
  Go to https://github.com/NOAA-GFDL/MOM6-examples/wiki and look at "Getting Started".

  Installation, compilation and running are platform specific operations for which we can only provide templates (as is
  done in on the wiki) but for which MOM6 developers cannot possibly support since every platform is different. Normally
  a user needs to know where libraries (such as netcdf and MPI) and compilers are on their system but once these have
  been established the documented compile process can be adapted to the local system.

User guide
  The `Read the Docs manual <https://mom6.readthedocs.io/en/main/>`_ provides
  a high-level overview of the model, its scientific formulation, and its use.

  The user guide is written in reStructuredText (.rst files) that reside in ``docs/`` of the `MOM6 source code <https://github.com/mom-ocean/MOM6>`_.
   The site renders these sources directly with Docutils, using Pybtex for citations and the
   bibliography. It requires neither Sphinx nor Doxygen.

   The longer scientific chapters are maintained as ``.rst`` files alongside
   the Fortran source under ``src/``. These replace the former ``.dox`` files
   and are included in the unified site. MathJax renders their equations
   in the browser; the site bundles it for offline use.

API reference
  The `experimental documentation site <https://marshallward.github.io/MOM6/>`_
  uses Flint to extract module and procedure documentation from Fortran source.
  Its API reference includes linked caller and callee diagrams, grouped by
  module, with expandable text views. These diagrams show possible calls,
  including private implementation procedures, rather than execution order.

  Read the Docs uses the same builder; older published versions may retain the
  former Doxygen-based API reference.

Repository policies
  Policies governing how the repositories are organized and operated live at https://github.com/NOAA-GFDL/MOM6-examples/wiki/MOM6-repository-policies.

Developer guide
  Beyond the API reference above, developer specific wiki pages are attached to the `MOM6 code repository <https://github.com/NOAA-GFDL/MOM6/wiki>`_.
