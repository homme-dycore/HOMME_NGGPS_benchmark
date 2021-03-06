#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([HOMME], [1.1.0], [jgarcia@ucar.edu],[homme], [http://www.homme.ucar.edu])
AC_CONFIG_AUX_DIR([../../../helpers]) 
AM_INIT_AUTOMAKE
#AC_CONFIG_HEADERS([config.h],[mv config.h config.hx && awk '/\/\*/,/\*\//{sub(/^/,"!")};{print}' config.hx >config.h && rm config.hx])
AC_CONFIG_HEADERS([config.h],[mv config.h config.hx && awk '/\/\*/,/\*\//{sub(/^\/\*.*/,"")};{print}' config.hx >config.h && rm config.hx])

# Set Defines that control the way the source units look after preprocessing.
AC_DEFINE([_MPI],[1],[Using MPI])
AC_DEFINE([_PRIM],[1],[Using PRIMITIVE])
AC_DEFINE([_PRESTART],[1],[Using _PRESTART])
AC_DEFINE([_WK_GRAD],[1],[_WK_GRAD])
AC_DEFINE([SPHEREW],[1],[SPHEREW])

# Checks for programs.
AC_PROG_RANLIB
AC_PROG_FC
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lblas':
#AC_CHECK_LIB([blas], [main])
# FIXME: Replace `main' with a function in `-llapack':
#AC_CHECK_LIB([lapack], [main])
# FIXME: Replace `main' with a function in `-lnetcdf':
#AC_CHECK_LIB([netcdf], [main])
# FIXME: Replace `main' with a function in `-lpio':
#AC_CHECK_LIB([pio], [main])
# FIXME: Replace `main' with a function in `-lpnetcdf':
#AC_CHECK_LIB([pnetcdf], [main])
# FIXME: Replace `main' with a function in `-ltiming':
#AC_CHECK_LIB([timing], [main])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.


# CHECK COMPILER NAME
if echo $FC | grep xlf >/dev/null 2>&1; then
       echo "Fortran Compiler is XLF"
       MY_DEFS=-WF,-DHAVE_CONFIG_H
   else
      MY_DEFS=-DHAVE_CONFIG_H  
fi

AC_ARG_ENABLE([energy-diagnostics],
            [AS_HELP_STRING([--enable-energy-diagnostics],
              [set the build so energy diagnostics are included])],
            [enable_energy_diagnostics=yes],
            [enable_energy_diagnostics=no])

AC_ARG_ENABLE([mesh-refine],
            [AS_HELP_STRING([--enable-mesh-refine],
              [set the build so Exodus library is included])],
            [enable_mesh_refine=yes],
            [enable_mesh_refine=no])

AC_ARG_ENABLE([bgl],
            [AS_HELP_STRING([--enable-bgl],
              [set the build so it is taylored for a BGL system])],
            [enable_bgl=yes],
            [enable_bgl=no])

AC_ARG_ENABLE([pio],
            [AS_HELP_STRING([--enable-pio],
              [set the build so PIO is used instead of default which is PIO_INTERP])],
            [enable_pio=yes],
            [enable_pio=no])

HOMME_LDFLAGS=""
HOMME_INCLUDE=""
HOMME_LIBS="-lnetcdf"

# defines if we want energy diagnostics in this build
AS_IF([test "x$enable_energy_diagnostics" = "xyes"], [
  dnl Do the stuff needed for enabling the feature
            echo 'Adding energy diagnostics in the config.h file'
            AC_DEFINE([ENERGY_DIAGNOSTICS],[1],[Energy Diagnostics])  
],[])

# defines if we want to include Exodus library in this build
AS_IF([test "x$enable_mesh_refine" = "xyes"], [
  dnl Do the stuff needed for enabling the feature
             echo 'Setting CPP macro to enable mesh refine'
            AC_DEFINE([MESH],[1],[using mesh refine]) 
],[])

# defines if we want to  configure this build for a BGL system
AS_IF([test "x$enable_bgl" = "xyes"], [
  dnl Do the stuff needed for enabling the feature
            echo 'Including _BGL in the config.h file'
            AC_DEFINE([_BGL],[1],[BGL system])  
],[])

# defines if we want pio in this build
AS_IF([test "x$enable_pio" = "xyes"], [
  dnl Do the stuff needed for enabling the feature
            echo 'Using PIO for input/output'
            AC_DEFINE([PIO],[1],[Basic PIO])  
],[
            echo 'Using PIO_INTERP for input/output'       
            AC_DEFINE([PIO_INTERP],[1],[PIO with interpolation]) 
])

AC_ARG_WITH(netcdf, [  --with-netcdf=PATH       specify prefix directory for installed NetCDF library.])

if test "x$with_netcdf" != x; then
  HOMME_LDFLAGS="-L$with_netcdf/lib $HOMME_LDFLAGS"
  HOMME_INCLUDE="-I$with_netcdf/include $HOMME_INCLUDE"
fi

AC_ARG_WITH(pnetcdf, [  --with-pnetcdf=PATH      specify prefix directory for installed pnetcdf library.])

if test "x$with_pnetcdf" != x; then
  HOMME_LDFLAGS="-L$with_pnetcdf/lib $HOMME_LDFLAGS"
  HOMME_LIBS="-lpnetcdf $HOMME_LIBS"
fi

AC_SUBST([HOMME_LDFLAGS])
AC_SUBST([HOMME_INCLUDE])
AC_SUBST([HOMME_LIBS])

AC_SUBST([MY_DEFS])
AC_SUBST([LOCAL_BLAS_LOCATION])
AC_SUBST([LOCAL_LAPACK_LOCATION])

AC_ARG_WITH(arch, [  --with-arch=<Linux, Darwin>  Specify you want to build for any of these specific architectures.])

if test "x$with_arch" != x; then
   if test "$with_arch" == "Linux"; then
      echo "Setting to use Linux architecture"   
      AC_DEFINE([_Linux],[1],[Linux System])
   else 
      if test "$with_arch" == "Darwin"; then
         echo "Setting to use Darwin architecture"   
         AC_DEFINE([_Darwin],[1],[Darwin System])
      else 
         echo "====>Unknown architecture $with_arch, ignoring directive"
      fi
   fi
fi

if test "${NP+set}" = set; then
  echo "variable NP is set to: $NP"
AC_SUBST([NP]) 
AC_DEFINE_UNQUOTED([NP],[${NP}],[Number of POINTS in the element (Default value is 4)])
else
  echo "variable NP is not set, using default value 4"
  AC_DEFINE_UNQUOTED([NP],[4],[Number of POINTS in the element (Default value is 4)])
fi 

if test "${NC+set}" = set; then
  echo "variable NC is set to: $NC"
AC_SUBST([NC]) 
AC_DEFINE_UNQUOTED([NC],[${NC}],[Number of cells NC*NC in each element (Default value NC=4])
else
  echo "variable NC is not set, using default value NC=4"
  AC_DEFINE_UNQUOTED([NC],[4],[Number of cells NC*NC in each element (Default value NC=4)])
fi

if test "${PLEV+set}" = set; then
  echo "variable PLEV is set to: $PLEV"
AC_SUBST([PLEV]) 
AC_DEFINE_UNQUOTED([PLEV],[${PLEV}],[Number of POINTS in the element (Default value is 26)])
else
  echo "variable PLEV is not set, using default value 26"
  AC_DEFINE_UNQUOTED([PLEV],[26],[Number of POINTS in the element (Default value is 26)])
fi 

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
