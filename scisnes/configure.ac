dnl SCISNES - Scilab Interface to the SNES solver
dnl Copyright (C) 2009 Yann COLLETTE
dnl
dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl (at your option) any later version.
dnl
dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

AC_PREREQ(2.63)
AC_INIT([scisnes], [1.0.0], [http://code.google.com/p/scilab-mip])

m4_include([m4/petsc.m4])
m4_include([m4/scilab.m4])

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR(sci_gateway/cpp/parameters.cpp)
AC_CONFIG_HEADER(sci_gateway/cpp/config.h)
AM_INIT_AUTOMAKE

dnl Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

if test x"$ac_cv_prog_cc_stdc" = x'no'; then
   AC_MSG_WARN([[found the $CC compiler but it is not ANSI-C compliant.]])
fi

AC_PROG_LIBTOOL

dnl Check to see if building shared libraries
libtool_build_shared_libs='no'
if ./libtool --config | grep 'build_libtool_libs=yes' > /dev/null
then
  libtool_build_shared_libs='yes'
fi

dnl check for Scilab
AC_CHECK_SCILAB

dnl Check for petsc
AC_CHECK_PETSC

dnl Checks for header files. 
AC_HEADER_STDC
if test x"$ac_cv_header_stdc" = x'no'; then
   AC_MSG_WARN([[this package uses ANSI-C headers but none were found.]])
fi
AC_CHECK_HEADERS([string.h])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_HEADER_STDBOOL

dnl Checks for library functions.
AC_CHECK_FUNCS([malloc])

AC_CONFIG_FILES([Makefile builder.sce])

AC_OUTPUT

AC_MSG_RESULT([
=====================================================
Configuration:
	CFLAGS               = $CFLAGS
	LDFLAGS              = $LDFLAGS
  	PETSC_CXXFLAGS       = $PETSC_CXXFLAGS
   	PETSC_LDADD          = $PETSC_LDADD
        USES_PETSC           = $USES_PETSC
	SCIEXE               = $SCIEXE
	SCILAB_VERSION_MAJOR = $SCILAB_VERSION_MAJOR
	SCILAB_VERSION_MINOR = $SCILAB_VERSION_MINOR
	SCILAB_VERSION_MICRO = $SCILAB_VERSION_MICRO
	PACKAGE_NAME         = $PACKAGE_NAME
	PACKAGE_VERSION      = $PACKAGE_VERSION
=====================================================
Please run 'make' to compile the toolbox. 
]) 
