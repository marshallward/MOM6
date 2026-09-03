#ifndef OMP_LOOP_BIND_COMPAT_H_
#define OMP_LOOP_BIND_COMPAT_H_

! This macro conditionally applies the bind(teams,parallel) clause of an
! OpenMP loop construct if supported by the compiler.
#ifdef HAVE_FC_OMP_LOOP_BIND
#define LOOP_BIND_TEAMS_PARALLEL bind(teams,parallel)
#else
#define LOOP_BIND_TEAMS_PARALLEL
#endif

#endif
