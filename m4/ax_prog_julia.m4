AC_DEFUN([AX_PROG_JULIA],[dnl
AC_ARG_VAR([JULIA],[Path to Julia executable.])
AC_PATH_PROG([JULIA], [julia], [no])
if test "$JULIA" = no; then
    AC_MSG_WARN([Julia not found. This is optional and only needed for running ./test scripts.
The path to a valid Julia executable may be specified with ./configure JULIA=path-to-julia])
fi;dnl
])
