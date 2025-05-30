dnl AX_FC_CHECK_MODULE(MODULE,
dnl                    [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND],
dnl                    [OTHER-FCFLAGS])
dnl
dnl This macro checks if a Fortran module is available to the compiler.
dnl
dnl The fourth argument (optional) allows for specification of supplemental
dnl FCFLAGS arguments.  This would primarily be used to test additional
dnl paths (typically using -I) for the module file.
dnl
dnl Results are cached in the ax_fc_cv_mod_MODULE variable.
dnl
AC_DEFUN([MOM6_FC_CHECK_MODULE],
[
  AC_LANG_ASSERT([Fortran])
  AS_VAR_PUSHDEF([mom6_fc_chk_mod_result], [mom6_cv_fc_mod_$1])
  AC_CACHE_CHECK([if $FC can use module $1], [mom6_cv_fc_mod_$1],[
    mom6_fc_chk_mod_save_FCFLAGS=$FCFLAGS
    test -n "$4" && FCFLAGS="$4 $FCFLAGS"
    AC_COMPILE_IFELSE(
        [AC_LANG_PROGRAM([],[use $1])],
        [AS_VAR_SET([mom6_fc_chk_mod_result], [yes])],
        [AS_VAR_SET([mom6_fc_chk_mod_result], [no])]
    )
    FCFLAGS=$mom6_fc_chk_mod_save_FCFLAGS
  ])
  AS_VAR_IF([mom6_fc_chk_mod_result], [yes], [$2], [$3])
  AS_VAR_POPDEF([mom6_fc_chk_mod_result])
])
