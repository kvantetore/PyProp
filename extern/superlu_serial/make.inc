############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
-include ../../../core/makefiles/Makefile.include
-include ../../../Makefile.platform
PLAT = _pyprop

#
#  The name of the libraries to be created/linked to
#
SUPERLULIB  = ../lib/libsuperlu_serial.a
BLASDEF	   	= -DUSE_VENDOR_BLAS
BLASLIB    	= $(LAPACK_LIBS)
LIBS		= $(SUPERLULIB) $(BLASLIB)

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib

CC           = $(MPICXX) -x c
CFLAGS       = $(CXXFLAGS) $(PIC)
FORTRAN	     = $(FC)
FFLAGS       = $(FORTRAN_FLAGS) $(PIC)
LOADER       = $(MPICXX)
LOADOPTS     = $(LDFLAGS)
NOOPTS       = $(PIC)
#
#  C preprocessor defs for compilation for the Fortran interface
#  (-DNoChange, -DAdd_, -DAdd__, or -DUpCase)
#
CDEFS        = -DAdd_
#
# The directory in which Matlab is installed
#


