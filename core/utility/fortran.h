#ifndef FORTRAN_H
#define FORTRAN_H

/*
 * Tools for calling fortran routines
 */

/* Fortran name mangler
 * On most platforms fortran names start with an underscore. However, this is (of course)
 * not the case on AIX
 */
#ifdef FORTRAN_NAME_NOUNDERSCORE
#define FORTRAN_NAME(x) x
#else
#define FORTRAN_NAME(x) x ## _
#endif

#endif

