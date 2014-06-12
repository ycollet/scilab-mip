
AC_PREREQ([2.52])

dnl @synopsis AX_TRILINOS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for libraries that implement the Trilinos ML
dnl and AztecOO packages.
dnl On success, it sets the TRILINOS_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with Trilinos, you should link with:
dnl
dnl     $TRILINOS_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order. LAPACK_LIBS is the output variable of the ACX_LAPACK
dnl macro, called automatically. BLAS_LIBS is the output variable of the
dnl ACX_BLAS macro, called automatically. FLIBS is the output variable
dnl of the AC_F77_LIBRARY_LDFLAGS macro (called if necessary by
dnl ACX_BLAS), and is sometimes necessary in order to link with F77
dnl libraries. Users will also need to use AC_F77_DUMMY_MAIN (see the
dnl autoconf manual), for the same reason.
dnl
dnl The user may also use --with-trilinos=DIR in order to specifiy
dnl the Trilinos install directory.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if the necessary
dnl libraries are found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if they are not found.
dnl
dnl @author Uche Mennel <umennel@student.ethz.ch>
dnl
AC_DEFUN([AX_TRILINOS], [
AC_REQUIRE([AX_LAPACK])
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_trilinos_ok=yes

acx_trilinos_save_LIBS="$LIBS"
TRILINOS_LIBS="-laztecoo -ltriutils -lteuchos -lepetraext -lepetra"
EXTRA_LIBS="$PARMETIS_LIBS $SUPERLU_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS $LIBS"


AC_ARG_WITH(trilinos,
AC_HELP_STRING([--with-trilinos],
[Specify Trilinos directory.  For example, --with-trilinos=/opt/trilinos]),
[
TRILINOS_INCLUDE="-I$withval/include"
TRILINOS_LIBS="-L$withval/lib $TRILINOS_LIBS"
CPPFLAGS="$CPPFLAGS -I$withval/include"
LDFLAGS="$LDFLAGS -L$withval/lib"
],)

# Find Trilinos libraries
AC_MSG_CHECKING(for AztecOO in Trilinos)
LIBS="$TRILINOS_LIBS $EXTRA_LIBS"
AC_TRY_LINK([
class AztecOO;
],
[AztecOO solver();],
[echo "yes"],
[acx_trilinos_ok=no; echo "no"])


AC_MSG_CHECKING(whether Ifpack was enabled in Trilinos)
LIBS="$TRILINOS_LIBS -lifpack $LIBS"
AC_TRY_LINK([
class Ifpack;
],
[Ifpack factory();],
[AC_DEFINE(HAVE_IFPACKLIB, [], [ifpack]) TRILINOS_LIBS="$TRILINOS_LIBS -lifpack"; echo "yes"],
[echo "no"])


AC_MSG_CHECKING(wether Amesos was enabled in Trilinos)
LIBS="$TRILINOS_LIBS -lamesos $EXTRA_LIBS"
AC_TRY_LINK([
class Amesos;
],
[Amesos factory();],
[AC_DEFINE(HAVE_AMESOSLIB, [], [amesos]) TRILINOS_LIBS="$TRILINOS_LIBS -lamesos"; echo "yes"],
[echo "no"])


AC_MSG_CHECKING(wether ML was enabled in Trilinos)
LIBS="$TRILINOS_LIBS -lml $EXTRA_LIBS"
AC_TRY_LINK([
namespace ML_Epetra {
class MultiLevelPreconditioner;
}
class Epetra_RowMatrix;
],
[Epetra_RowMatrix *A;
ML_Epetra::MultiLevelPreconditioner Prec();],
[AC_DEFINE(HAVE_MLLIB, [], [ml]) TRILINOS_LIBS="$TRILINOS_LIBS -lml"; echo "yes"],
[echo "no"])

LIBS="$acx_trilinos_save_LIBS"

AC_SUBST(TRILINOS_LIBS)
AC_SUBST(TRILINOS_INCLUDE)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_trilinos_ok" = xyes; then
        $1
        :
else
        acx_trilinos_ok=no
        $2
fi
])
 dnl ACX_TRILINOS
