AC_DEFUN([AC_CHECK_PETSC],
[
  AH_TEMPLATE(HAVE_PETSC,
   [Defined to 1 if PETSc is present on the system])

  AC_ARG_WITH(petsc_prefix,
		AC_HELP_STRING([--with-petsc-prefix=DIR],[Set the path to the PETSC]),
		[with_petsc_prefix=$withval],
		[with_petsc_prefix='yes']
		)
 
  REQUIRED_PETSC_MAJOR=3
  REQUIRED_PETSC_MINOR=0
  PETSC_VERSION_MAJOR=0
  PETSC_VERSION_MINOR=0

  REQUIRED_PETSC=`echo "$REQUIRED_PETSC_MAJOR.$REQUIRED_PETSC_MINOR"`

  petscinclude=""
  petsclibpath=""

  dnl look for petsc in standard pathes
  
  if test "x$with_petsc_prefix" != "xyes"
  then
    path_list="/usr/ /usr/local/ $with_petsc_prefix/"
  else
    path_list="/usr/ /usr/local/"
  fi
 
  for i in $path_list
  do 
    if test -e "$i/include/petscversion.h"
    then
      petscinclude="$i/include"

      if test "x`ls $i/lib/libpetsc.* 2> /dev/null`" != "x"
      then
        petsclibpath="$i/lib"
      fi
    fi
  done

  if test "x$petscinclude" != "x"
  then
    PETSC_CXXFLAGS="-I$petscinclude"
  fi
  if test "x$petsclibpath" != "x"
  then
    PETSC_LDADD="-L$petsclibpath -lpetsccontrib -lpetscdm -lpetscksp -lpetscmat -lpetscsnes -lpetscvec -lpetscts -lpetsc"
  fi

  dnl Compare MINIMUM-VERSION with libparted version
  config_major_version=3
  config_minor_version=0

  
  AC_TRY_RUN([
  #include <stdio.h>
  #include "$petscinclude/petscversion.h"

  int main ()
  {
    int major = PETSC_VERSION_MAJOR;
    int minor = PETSC_VERSION_MINOR;

    if ((major > $config_major_version) ||
       ((major == $config_major_version) && (minor >= $config_minor_version))) {
         return 0;
       } else {
         printf(" *** Looked for petsc version >= %d.%d, found %d.%d\n",
                $config_major_version, $config_minor_version, major, minor);
         return 1;
       }
  }
  ],
  AC_MSG_RESULT(yes); USES_PETSC="yes" , 
  AC_MSG_RESULT(no) ; USES_PETSC="no"; $3,
  [echo $ac_n "cross compiling; assumed OK... $ac_c"])

  dnl AC_DEFINE(HAVE_PETSC)

  AC_SUBST(PETSC_CXXFLAGS)
  AC_SUBST(PETSC_LDADD)
  AC_SUBST(USES_PETSC)
])
