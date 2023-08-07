AC_INIT([Kestrel], [1.0], [Mark.Woodhouse@bristol.ac.uk, J.Langham@bristol.ac.uk])
AC_CONFIG_SRCDIR([src/TimeStepper.f90])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Default install is local (in src/../)
AC_PREFIX_DEFAULT([..])

# Checks for compilers
AC_PROG_FC
AC_FC_LIBRARY_LDFLAGS
AC_PROG_CXX

# Check c++17, which we require.
AX_CXX_COMPILE_STDCXX([17], [noext], [mandatory])

# Need to find either the filesystem header, or boost
AC_LANG(C++)
AC_CHECK_HEADERS([filesystem], [], 
    # Failed to find filesystem, so check for boost
    [
     AX_BOOST_BASE([1.46.0], [], 
       [
        AC_MSG_ERROR([Either the C++17 filesystem library
        is required (e.g. gcc --version >= 9), or we need Boost.])])
     AX_BOOST_SYSTEM
     AX_BOOST_FILESYSTEM
     ]
)
AM_CONDITIONAL([WANT_BOOST], [test x$want_boost = xyes])

AX_LIB_GDAL([2.2.0])

AX_LIB_PROJ()

AX_LIB_NETCDF4([serial])

AX_PROG_JULIA()
AM_CONDITIONAL([HAVE_JULIA], [test x$JULIA != x])

AC_ARG_ENABLE([debug],
[  --enable-debug          turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-debug]) ;;
esac],[debug=false])
AM_CONDITIONAL([DEBUG], [test x$debug = xtrue])

AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT
