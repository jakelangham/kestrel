AC_PREREQ([2.68])
m4_define([DEFAULT_VERSION], [v1.1.1])
AC_INIT([Kestrel],
        [m4_esyscmd_s([git --version > /dev/null 2>&1 && git describe --tags || echo DEFAULT_VERSION])],
        [Mark.Woodhouse@bristol.ac.uk, J.Langham@bristol.ac.uk])
AC_CONFIG_SRCDIR([src/TimeStepper.f90])
AC_CONFIG_MACRO_DIRS([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Default install is local (in src/../)
AC_PREFIX_DEFAULT([..])

# Checks for compilers
AC_PROG_FC
AC_FC_LIBRARY_LDFLAGS
AC_PROG_CXX

# Detect builds on linux
AC_CANONICAL_HOST
AM_CONDITIONAL([LINUX], [test x$host_os = xlinux-gnu])

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

AC_ARG_ENABLE([parallel],
[  --enable-parallel          turn on parallel loops],
[case "${enableval}" in
  yes) parallel=true ;;
  no)  parallel=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-parallel]) ;;
esac],[parallel=false])
AM_CONDITIONAL([PARAL], [test x$parallel = xtrue])

if test x$parallel = xtrue; then
#   FC_PARALLEL_FLAGS="-ftree-parallelize-loops=$threads -fopenmp"
  FC_PARALLEL_FLAGS="-fopenmp"
else
  FC_PARALLEL_FLAGS=""
fi

AC_SUBST([FC_PARALLEL_FLAGS])

AC_ARG_ENABLE([opt-report],
[  --enable-opt-report        Show missed optimization and vectorization reports],
[case "${enableval}" in
  yes) opt_report=true ;;
  no)  opt_report=false ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-opt-report]) ;;
esac],
[opt_report=false])

# Set Fortran flags accordingly
if test x$opt_report = xtrue; then
  FC_OPT_REPORT_FLAGS="-fopt-info-missed -fopt-info-vec-missed"
else
  FC_OPT_REPORT_FLAGS=""
fi

AC_SUBST([FC_OPT_REPORT_FLAGS])

AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT
