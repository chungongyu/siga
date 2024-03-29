#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([siga], [0.10], [siga@ict.ac.cn], [siga], [http://bioinfo.ict.ac.cn/siga])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_AR
AC_PROG_LIBTOOL 
AC_OPENMP
AC_CHECK_PROG([with_help2man], [help2man], [yes], [no])
AM_CONDITIONAL([HELP2MAN], [test "x${with_help2man}" = "xyes"])
AX_BOOST_BASE([1.41], [], AC_MSG_ERROR([Could not find a useful version of boost]))
AX_BOOST_FILESYSTEM
AX_BOOST_IOSTREAMS
AX_BOOST_REGEX
AX_BOOST_SYSTEM
AX_BOOST_UNIT_TEST_FRAMEWORK

AX_CXX_COMPILE_STDCXX_11([noext],[mandatory])

PKG_CHECK_MODULES([RAPIDJSON], [RapidJSON >= 1.0.0])
PKG_CHECK_MODULES([LOG4CXX], [liblog4cxx >= 0.10.0])
AC_ARG_WITH([tcmalloc], AS_HELP_STRING([--with-tcmalloc], [Build with the `TCMalloc` (default=no)]))
AS_IF([test "x${with_tcmalloc}" = "xyes"], [
    PKG_CHECK_MODULES([TCMALLOC], [libtcmalloc >= 2.1])
])
PKG_CHECK_MODULES([MLPACK], [mlpack >= 2.3.1], [
    AC_DEFINE([HAVE_MLPACK], [1], [have mlpack])], [
])

AC_ARG_ENABLE([gprof], [AS_HELP_STRING([--enable-gprof], [enable internal support for `gprof` (default=no)])],
[case "${enableval}" in
    yes)
        CFLAGS="$CFLAGS -pg"
        CXXFLAGS="$CXXFLAGS -pg"
        ;;
    no) ;;
    *) AC_MSG_ERROR(bad value ${enableval} for --enable-gprof);;
esac])

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h omp.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset pow sqrt])

AC_CONFIG_FILES([Makefile
                 3rdparty/Makefile
                 benchmark/Makefile
                 docs/Makefile
                 docs/_layouts/Makefile
                 examples/Makefile
                 m4/Makefile
                 src/Makefile
                 test/Makefile])
AC_CONFIG_LINKS([docs/Gemfile:docs/Gemfile])
AC_CONFIG_SUBDIRS()
AC_OUTPUT
