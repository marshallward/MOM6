dnl MOM6_FC_OPENMP_DTYPE_ALLOCS
dnl
dnl Test OpenMP mapping of derived types with allocatable components.
dnl Call with OpenMP enabled in FCFLAGS.
dnl
AC_DEFUN([MOM6_FC_OPENMP_DTYPE_ALLOCS], [
  AC_REQUIRE([AC_PROG_FC])
  AC_LANG_PUSH([Fortran])
  AC_CACHE_CHECK(
    [whether $FC supports OpenMP mapping of derived types with allocatable components],
    [mom6_cv_fc_openmp_dtype_allocs], [
      AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([], [
dnl ---
      implicit none
      type :: container
        real, allocatable :: field(:)
      end type container
      type(container) :: obj

      allocate(obj%field(4))
      obj%field(:) = 1.

      !$omp target enter data map(to: obj)
      !$omp target enter data map(to: obj%field)
      !$omp target exit data map(delete: obj%field)
      !$omp target exit data map(delete: obj)

      deallocate(obj%field)
dnl ---
        ])
      ], [
        mom6_cv_fc_openmp_dtype_allocs=yes
      ], [
        mom6_cv_fc_openmp_dtype_allocs=no
      ])
    ])
  AS_IF([test "$mom6_cv_fc_openmp_dtype_allocs" = yes], [
    AC_DEFINE([HAVE_FC_OPENMP_DTYPE_ALLOCS], [1],
      [Define if OpenMP supports mapping derived types with allocatable components.])
  ])
  AC_LANG_POP([Fortran])
])
