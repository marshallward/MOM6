dnl Determine the flag required to force 64-bit reals.
dnl
dnl Many applications do not specify the kind of its real variables, even
dnl though the code may intrinsically require double-precision.  Most compilers
dnl will also default to using single-precision (32-bit) reals.
dnl
dnl This test determines the flag required to set reals without explcit kind to
dnl 64-bit double precision floats.  Ideally, we also desire to leave any
dnl `DOUBLE PRECISION` variable as 64-bit.  But this does not appear to always
dnl be possible, such as in NAG Fortran (see below).
dnl
dnl This does not test if the behavior of integers is changed; for example,
dnl Cray's Fortran wrapper's -default will double both.  This is addressed by
dnl avoiding any flags with affect integers, but this should still be used with
dnl some care.
dnl
dnl   GCC               -fdefault-real-8, -fdefault-double-8
dnl   AMD (flang)       -fdefault-real-8
dnl   [Common alias]    -r8
dnl   Intel             -real-kind 64
dnl   PGI/Nvidia        -Mr8
dnl   Cray              -s real64
dnl   NAG               -double
dnl
dnl
dnl NOTE:
dnl   - Many compilers accept -r8 for real and double precision sizes, but
dnl     several compiler-specific options are also provided.
dnl
dnl   - -r8 in NAG will attempt to also set double precision to 16 bytes if
dnl     available, which is generally undesired.
dnl
dnl     Additionally, the -double flag, which doubles *all* types, appears to
dnl     be the preferred flag here.
dnl
dnl     Neither flag describes what we actually want, but we include it here
dnl     as a last resort.
dnl
AC_DEFUN([MOM6_FC_REAL64], [
  AC_ARG_ENABLE([real64], [
    AS_HELP_STRING([--disable-real-64], [do not use 64-bit reals])
  ])
  AC_ARG_ENABLE([strict-real64], [
    AS_HELP_STRING([--enable-strict-real64],
      [abort if no 8-byte reals are unsupported]
    )
  ])
  AC_ARG_VAR([REAL64_FCFLAGS], [User-defined Fortran flag to set default reals to 64-bit])
  AS_IF([test "x$enable_real64" != xno], [
    AS_IF([test -n "$REAL64_FCFLAGS"], [
      AC_MSG_NOTICE([Using user-defined REAL64_FCFLAGS: $REAL64_FCFLAGS])
    ], [
      AC_CACHE_CHECK([for $FC option to force 64-bit reals],
        [mom6_cv_prog_fc_real64], [
          mom6_cv_prog_fc_real64="unsupported"
          ac_fc_real64_FCFLAGS_save=${FCFLAGS}
          for ac_flag in none \
            -fdefault-real-8 \
            "-fdefault-real-8 -fdefault-double-8" \
            -r8 \
            "-real-kind 64" \
            -Mr8 \
            "-s real64" \
            -double
          do
            test "$ac_flag" != none \
              && FCFLAGS="$ac_fc_real64_FCFLAGS_save $ac_flag"
            AC_LINK_IFELSE([
              AC_LANG_PROGRAM([], [
dnl ---
     real :: x(4)
     double precision :: y(4)
     integer, parameter :: &
       m = merge(1, 0, kind(x(1)) == selected_real_kind(15, 307)), &
       n = merge(1, 0, kind(y(1)) == selected_real_kind(15, 307))
     print *, x(::m)
     print *, y(::n)
dnl ---
              ])
            ], [
              mom6_cv_prog_fc_real64="$ac_flag"
              break
            ])
          done
          FCFLAGS=$ac_fc_real64_FCFLAGS_save
        ]
      )
      AS_CASE([$mom6_cv_prog_fc_real64],
        [none], [
          mom6_cv_prog_fc_real64="none needed"
          REAL64_FCFLAGS=""
        ],
        [unsupported], [
          AS_IF([test "x$enable_strict_real64" = xyes], [
            AC_MSG_ERROR(
              [No known flag found to force 64-bit reals; aborting]
            )
          ], [
            AC_MSG_WARN(
              [No known flag found to force 64-bit reals; using defaults]
            )
            REAL64_FCFLAGS=""
          ])
        ],
        [REAL64_FCFLAGS=$mom6_cv_prog_fc_real64]
      )
    ])
  ])
  AC_SUBST([REAL64_FCFLAGS])
])
