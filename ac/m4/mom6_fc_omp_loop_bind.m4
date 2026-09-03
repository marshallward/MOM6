dnl MOM6_FC_OMP_LOOP_BIND
dnl
dnl Determine if the compiler supports the `bind` clause of the OpenMP `loop`
dnl construct, e.g. `!$omp loop bind(teams,parallel)`.  This clause was
dnl introduced in OpenMP 5.1 and is not yet universally supported.
dnl
dnl Results are cached in the `mom6_cv_fc_omp_loop_bind` variable.
dnl
AC_DEFUN([MOM6_FC_OMP_LOOP_BIND], [
  AC_REQUIRE([AC_PROG_FC])
  AC_LANG_PUSH([Fortran])
  AC_CACHE_CHECK([whether $FC supports omp loop bind(teams,parallel)],
    [mom6_cv_fc_omp_loop_bind], [
      mom6_cv_fc_omp_loop_bind=no
      mom6_fc_omp_loop_bind_save_FCFLAGS=$FCFLAGS
      FCFLAGS="$OPENMP_FCFLAGS $FCFLAGS"
      AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([], [
dnl ---
      !\$omp target teams
      !\$omp loop bind(teams,parallel)
      do i=1,2
      end do
      !\$omp end target teams
dnl ---
        ])
      ], [
        mom6_cv_fc_omp_loop_bind=yes
      ])
      FCFLAGS=$mom6_fc_omp_loop_bind_save_FCFLAGS
    ]
  )
  AS_IF([test "x$mom6_cv_fc_omp_loop_bind" = "xyes"], [
    AC_DEFINE([HAVE_FC_OMP_LOOP_BIND], [1],
      [Define if the OpenMP loop construct supports the bind clause.])
  ])
  AC_LANG_POP([Fortran])
])
