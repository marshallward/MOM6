! The UNUSED macro generates a null operation which forces most compilers to
! recognize a variable as used, while also removing the statement during.
! optimization.
#ifdef ENABLE_UNUSED
#define UNUSED(x) associate(ignored_variable_ => x); end associate
#else
#define UNUSED(x)
#endif
