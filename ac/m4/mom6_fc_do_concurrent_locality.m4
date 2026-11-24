dnl MOM6_FC_DO_CONCURRENT_LOCALITY
dnl
dnl Determine if the compiler supports do-concurrent locality specifiers.
dnl LOCAL, LOCAL_INIT, and REDUCE are tested, but this should also include
dnl SHARED and DEFAULT(NONE).
dnl
dnl Results are cached in the `mom6_cv_fc_do_concurrent_locality` variable.
dnl
AC_DEFUN([MOM6_FC_DO_CONCURRENT_LOCALITY], [
  AC_REQUIRE([AC_PROG_FC])
  AC_LANG_PUSH([Fortran])
  AC_CACHE_CHECK([whether $FC supports do concurrent locality specifiers],
    [mom6_cv_fc_do_concurrent_locality], [
      mom6_cv_fc_do_concurrent_locality=no
      AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM([], [
dnl ---
      implicit none
      integer :: i, a, b, total, bins(2), largest
      logical :: found
      total = 0
      bins(:) = 0
      largest = 0
      found = .false.
      b = 2
      do concurrent(i=1:2) local(a) local_init(b) reduce(+:total,bins)
        a = i
        b = b * i
        total = total + a
        bins(:) = bins(:) + b
      end do
      do concurrent(i=1:2) reduce(max:largest) reduce(.or.:found)
        largest = max(largest, i)
        found = found .or. (i == 2)
      end do
dnl ---
        ])
      ], [
        mom6_cv_fc_do_concurrent_locality=yes
      ])
    ]
  )
  AS_IF([test "$mom6_cv_fc_do_concurrent_locality" = "yes"], [
    AC_DEFINE([HAVE_FC_DO_CONCURRENT_LOCALITY], [1],
      [Define if do concurrent supports locality specifiers.])
  ])
  AC_LANG_POP([Fortran])
])
