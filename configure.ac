AC_PREREQ([2.68])
m4_define([DEFAULT_VERSION], [v1.1.2])
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

# Check for -march=native -mtune=native,
# otherwise try -mcpu=native -mtune=native, or else give up. 
# N.B. These flags seem to conflict for some compilers so must be tested
# together.
AC_LANG(C++)
AX_CHECK_COMPILE_FLAG(-march=native -mtune=native,
                      [AC_SUBST(CXX_ARCH_FLAGS, "-march=native -mtune=native")],
                      [AX_CHECK_COMPILE_FLAG(-mcpu=native -mtune=native,
                        [AC_SUBST(CXX_ARCH_FLAGS, "-mcpu=native -mtune=native")],
                       )])
AC_LANG(Fortran)
AX_CHECK_COMPILE_FLAG(-march=native -mtune=native,
                      [AC_SUBST(FC_ARCH_FLAGS, "-march=native -mtune=native")],
                      [AX_CHECK_COMPILE_FLAG(-mcpu=native -mtune=native,
                        [AC_SUBST(FC_ARCH_FLAGS, "-mcpu=native -mtune=native")],
                       )])

# OpenMP (multithreading) is off by default; turn on with --enable-openmp.
# OPENMP_FCFLAGS is used to compile the Fortran sources; OPENMP_CXXFLAGS is
# needed at link time, since we link with the C++ compiler. (AC_OPENMP only
# looks for these flags when enable_openmp != no.)
AC_ARG_ENABLE([openmp],
[  --enable-openmp         build with OpenMP multithreading (default: no)],
[case "${enableval}" in
  yes) enable_openmp=yes ;;
  no)  enable_openmp=no ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-openmp]) ;;
esac],[enable_openmp=no])
# AC_OPENMP declares its own (inert) --disable-openmp option, which would
# duplicate the one above in --help, so we suppress it.
m4_pushdef([AC_ARG_ENABLE], [])
AC_LANG(Fortran)
AC_OPENMP
AC_LANG(C++)
AC_OPENMP
m4_popdef([AC_ARG_ENABLE])
AS_IF([test "x$enable_openmp" = xyes],
      [AS_IF([test "x$ac_cv_prog_fc_openmp" = xunsupported || \
              test "x$ac_cv_prog_cxx_openmp" = xunsupported],
             [AC_MSG_ERROR([--enable-openmp was given, but the compilers do not support OpenMP])])
       AC_MSG_NOTICE([building with OpenMP ($OPENMP_FCFLAGS)])],
      [AC_MSG_NOTICE([building without OpenMP (use --enable-openmp to turn on)])])

AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT
