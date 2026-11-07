Source documentation completeness
=================================

``docs-and-style.yml`` checks whitespace and source documentation with
``flint check-docs``, then tests the builder and generates/validates the published site.
These are separate checks: an undocumented variable can still be rendered, but
is a source-policy failure. Neither check rewrites Fortran source.

Local command, from the repository root::

    flint check-docs src config_src/drivers/solo_driver \
        config_src/memory/dynamic_symmetric config_src/external \
        config_src/drivers/FMS_cap \
        -I config_src/memory/dynamic_symmetric \
        --defer-include src/framework/version_variable.h \
        --exclude src/equation_of_state/TEOS10 \
        --exclude pkg/GSW-Fortran --exclude pkg/CVMix-src

Install the Flint revision from ``requirements-api.txt`` first. The published
pin includes ``check-docs``, the declaration-discovery improvements, and the
scope-comment ownership fix: a member's trailing comment must not satisfy its
enclosing or following scope's description.

Scope
-----

The roots match the former Doxygen check, including both the solo driver and
FMS cap. They intentionally differ from site generation: files are checked
individually, so alternative drivers need not share a module namespace.
TEOS10 and dependency implementations are excluded, including symlink targets.

The default check requires descriptions for parsed concrete types,
named interfaces, module variables, and type components. Private entities are included, as in
the old ``EXTRACT_PRIVATE = YES`` configuration. Ordinary procedure-local
variables are excluded. Explicit ``!>``, ``!<``, and ``!!`` comments and
declaration groups qualify; empty markers and ordinary comments do not.
Descriptions shared by a declaration on one line count for all its entities;
a comment on a continued declaration does not automatically cover other lines.

Module/submodule, program, subroutine, function, and abstract-type descriptions
are also required by the intended
policy, but enforcement is deferred while the existing omissions are addressed.
Add ``--strict`` to include those descriptions. CI currently uses relaxed mode.
This does not skip their contents: variables, types, and interfaces inside
undocumented modules and procedures remain checked. Procedure descriptions are
strict-only regardless of public/private accessibility; the requirement is
deferred, not abandoned.
Abstract-type components remain checked in relaxed mode, as do concrete types
that extend an abstract base. Regression tests verify that ``--strict`` restores
each deferred description category. There are no permanent exclusions for them.

For initial adoption, CI also defers descriptions originating in the specific
header ``src/framework/version_variable.h`` using ``--defer-include``. The
declaration of ``version`` remains subject to the intended documentation policy;
``--strict`` restores this check too. Other headers and local declarations with
macro initializers are not exempt. Included declarations are reported using
their original header path and line, rather than the next line of the caller.

Private module parameters are also deferred to ``--strict``. Both properties
must hold: public parameters and private non-parameter variables are still
checked. The checker honors module default access, named ``public``/``private``
statements, and access attributes on the declaration.

The former ``WARN_NO_PARAMDOC = NO`` setting did not mandate documentation for
every dummy argument. Accordingly, this CI check does not enable Flint's optional
``--arguments`` mode. Specific type-bound procedure and finalization-procedure
descriptions are strict-only, like other procedures. Generic binding descriptions
are checked in relaxed mode, like named interfaces. Function-result descriptions,
units syntax, and prose quality are not checked. This policy is
explicit; it does not claim exact equivalence with Doxygen's heuristics or
coverage of unsupported parser constructs.

Dummy procedures declared with ``procedure(interface)`` or ``procedure()`` are
included in argument checking, with or without ``pointer``. Their declaration
comments and locations are retained. Module-level non-pointer procedure
declarations use the strict-only procedure-description policy; procedure pointers
stored in modules remain subject to the ordinary variable policy.

Diagnostics and adoption
------------------------

Missing descriptions are reported as ``path:line: missing-documentation: kind
qualified_name``. The command exits 1 for omissions and 2 for invalid inputs or
reported incomplete parses, including missing include files. Other parser
failures also terminate unsuccessfully. There is no suppressed-finding baseline
and no success status for an empty source selection.

The initial strict MOM6 run on 2026-11-07 checked 268 files and reported 231 candidates:
18 module variables, 183 subroutines, 21 functions, seven modules, one program,
and one derived type. Relaxed mode defers the seven modules and 183 subroutines,
leaving 41 initial findings. A subsequent macro-comment association fix removed
the false positive for ``posix::SIGUSR1``, leaving 40 candidates (17 module
variables, 21 functions, one program, and one type). Deferring program descriptions
to strict mode reduces the relaxed report to 39 candidates. Deferring private
and internal function descriptions removes another 18, leaving 21 candidates:
17 variables, three public functions, and one type. Deferring all function
descriptions leaves 18 candidates: 17 variables and one type. Deferring
abstract-type descriptions leaves 17 variable candidates. The check correctly exits
nonzero on this tree. Deferring ``version_variable.h`` removes eight repeated
findings, leaving nine module variables in ``MOM_file_parser_tests``. All nine
are private parameters; deferring that category leaves no relaxed-mode findings.
The strict report still exposes the backlog, and tests verify that its checks
remain active. Passing relaxed mode does not establish complete documentation.

The separate `symbol-inventory audit <symbol-inventory-comparison.rst>`_ compares
discovered declarations before policy exclusions. Its initial run fails the
Doxygen-subset requirement: procedure-pointer declarations and type-bound
bindings were among the gaps. Subsequent Flint procedure-pointer support recovers
the four module/component pointer declarations and checks their descriptions.
Type-bound binding support subsequently recovered 398 local binding declarations;
their descriptions are now subject to the procedure/generic policy above.
The completed entity-coverage audit found no remaining unaccounted-for Fortran
declarations in the selected sources and preprocessing configuration. Structural
representation differences, macro documentation, and comment-created namespaces
were reviewed and accepted. The raw comparison still fails pending normalization;
its allowlist remains empty and it is not a normal CI gate.

Adoption validation
-------------------

A clean installation of published Flint revision
``eb3d4187fd3398d314239f9d4c4faac738bc542a`` passed the relaxed check on all
268 selected files and all 18 MOM6 documentation/tooling tests. Full site
generation produced 340 HTML pages and 4,033 graphs, with zero unresolved
references and 161,589 checked links (126,940 cross-file).

Negative CLI checks exposed a scope-comment ownership bug in that revision.
The correction, published as ``2a789c8d5d62af62e77938a8de0f998861fffb32``
and now pinned in MOM6, passes Flint's full 329-test suite (one skipped), preserves
zero relaxed MOM6 findings, and reports 241 strict findings. CLI smoke checks
verify exits 0 for documented inputs, 1 for missing descriptions, and 2 for
missing includes; the CI ``pipefail``/``tee`` pipeline preserves failure status.
The strict count includes previously hidden scope-description omissions; the
adoption policy and its deferrals are unchanged.

Installation of the corrected published pin was verified independently of the
local Flint checkout. It passes the 268-file relaxed check, all 18 MOM6
documentation/tooling tests, and the strict scope-comment regression smoke check.

These validations used Python 3.14 locally. Hosted GitHub Actions and Read the
Docs runs using Python 3.13 remain to be verified.
