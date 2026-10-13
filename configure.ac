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
# OPENMP_FCFLAGS is used to compile the Fortran sources. Since we link with the
# C++ compiler, it also needs to link the Fortran compiler's OpenMP runtime:
# OPENMP_LINKFLAGS and OPENMP_LIBS. If the C++ compiler is GCC, this is just
# OPENMP_CXXFLAGS. Otherwise (e.g. Apple clang, which is needed to match
# Homebrew's GDAL on macOS, with Homebrew's gfortran), we link the GNU OpenMP
# runtime (libgomp) explicitly. (AC_OPENMP only looks for flags when
# enable_openmp != no.)
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
OPENMP_LINKFLAGS=
OPENMP_LIBS=
AS_IF([test "x$enable_openmp" = xyes],
      [AS_IF([test "x$ac_cv_prog_fc_openmp" = xunsupported],
             [AC_MSG_ERROR([--enable-openmp was given, but the Fortran compiler does not support OpenMP])])
       AS_IF([test "x$ac_cv_prog_cxx_openmp" != xunsupported && \
              ! $CXX --version 2>/dev/null | grep -qi clang],
             [OPENMP_LINKFLAGS=$OPENMP_CXXFLAGS],
             [AC_MSG_CHECKING([for the Fortran compiler's OpenMP runtime])
              gomp_dir=
              for gomp_lib in libgomp.dylib libgomp.so libgomp.a; do
                 gomp_path=`$FC -print-file-name=$gomp_lib 2>/dev/null`
                 AS_IF([test -f "$gomp_path"],
                       [gomp_dir=`AS_DIRNAME(["$gomp_path"])`; break])
              done
              AS_IF([test -z "$gomp_dir"],
                    [AC_MSG_RESULT([not found])
                     AC_MSG_ERROR([--enable-openmp was given, but could not find libgomp for $FC])])
              AC_MSG_RESULT([$gomp_dir])
              OPENMP_LINKFLAGS="-L$gomp_dir -Wl,-rpath,$gomp_dir"
              OPENMP_LIBS="-lgomp"
              save_LDFLAGS=$LDFLAGS; save_LIBS=$LIBS
              LDFLAGS="$LDFLAGS $OPENMP_LINKFLAGS"; LIBS="$OPENMP_LIBS $LIBS"
              AC_MSG_CHECKING([whether $CXX can link libgomp])
              AC_LINK_IFELSE([AC_LANG_PROGRAM([[extern "C" int omp_get_max_threads(void);]],
                                              [[return omp_get_max_threads() < 1;]])],
                             [AC_MSG_RESULT([yes])],
                             [AC_MSG_RESULT([no])
                              AC_MSG_ERROR([--enable-openmp was given, but $CXX cannot link $FC's OpenMP runtime (libgomp)])])
              LDFLAGS=$save_LDFLAGS; LIBS=$save_LIBS])
       AC_MSG_NOTICE([building with OpenMP ($OPENMP_FCFLAGS; linking with $OPENMP_LINKFLAGS $OPENMP_LIBS)])],
      [AC_MSG_NOTICE([building without OpenMP (use --enable-openmp to turn on)])])
AC_SUBST([OPENMP_LINKFLAGS])
AC_SUBST([OPENMP_LIBS])

# Link-time optimisation is off by default; turn on with --enable-lto. It lets
# the compiler inline small, frequently called routines across source files,
# which speeds up the solver by several percent. The Fortran and C++ object
# files are optimised together at link time, so this needs the Fortran and
# C++ compilers to be the same version of GCC (e.g. not Homebrew gfortran with
# Apple clang).
AC_ARG_ENABLE([lto],
[  --enable-lto            build with link-time optimisation (default: no)],
[case "${enableval}" in
  yes) enable_lto=yes ;;
  no)  enable_lto=no ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-lto]) ;;
esac],[enable_lto=no])
LTO_FLAGS=
AS_IF([test "x$enable_lto" = xyes],
      [AC_CACHE_CHECK([whether $FC and $CXX are the same version of GCC],
                      [kestrel_cv_lto_compatible],
                      [kestrel_cv_lto_compatible=no
                       fc_version=`$FC -dumpfullversion 2>/dev/null`
                       cxx_version=`$CXX -dumpfullversion 2>/dev/null`
                       AS_IF([test -n "$fc_version" && test "x$fc_version" = "x$cxx_version" && \
                              ! $CXX --version 2>/dev/null | grep -qi clang],
                             [kestrel_cv_lto_compatible=yes])])
       AS_IF([test "x$kestrel_cv_lto_compatible" != xyes],
             [AC_MSG_ERROR([--enable-lto needs the Fortran and C++ compilers to be the same version of GCC ($FC: ${fc_version:-unknown}, $CXX: ${cxx_version:-unknown})])])
       AC_LANG_PUSH([Fortran])
       AX_CHECK_COMPILE_FLAG([-flto=auto], [LTO_FLAGS="-flto=auto"],
                             [AX_CHECK_COMPILE_FLAG([-flto], [LTO_FLAGS="-flto"])])
       AC_LANG_POP([Fortran])
       AS_IF([test -z "$LTO_FLAGS"],
             [AC_MSG_ERROR([--enable-lto was given, but $FC does not support -flto])])
       AC_LANG_PUSH([C++])
       save_CXXFLAGS=$CXXFLAGS; save_LDFLAGS=$LDFLAGS
       CXXFLAGS="$CXXFLAGS $LTO_FLAGS"; LDFLAGS="$LDFLAGS $LTO_FLAGS"
       AC_MSG_CHECKING([whether $CXX can link with $LTO_FLAGS])
       AC_LINK_IFELSE([AC_LANG_PROGRAM([], [[return 0;]])],
                      [AC_MSG_RESULT([yes])],
                      [AC_MSG_RESULT([no])
                       AC_MSG_ERROR([--enable-lto was given, but $CXX cannot link with $LTO_FLAGS])])
       CXXFLAGS=$save_CXXFLAGS; LDFLAGS=$save_LDFLAGS
       AC_LANG_POP([C++])
       AC_MSG_NOTICE([building with link-time optimisation ($LTO_FLAGS)])])
AC_SUBST([LTO_FLAGS])

AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT
