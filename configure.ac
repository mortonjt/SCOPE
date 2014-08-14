#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([SCOPE], [1.0], [mortonjt@miamioh.edu])
AC_CONFIG_SRCDIR([src/SCOPE.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([no-define])
#m4_include([m4/ax_pthread.m4])
# Checks for programs.
AC_PROG_CXXCPP
AC_PROG_CXX
#AC_PROG_LIBTOOL

AC_LANG([C++])
AX_CXX_COMPILE_STDCXX_11
#AX_CHECK_COMPILE_FLAG([-Wall -Wextra -Wno-unused-variable -Wno-unused-but-set-variable \
#    -Wno-unused-parameter],\
#    [CXXFLAGS="$CXXFLAGS -Wall -Wextra -Wno-unused-variable -Wno-unused-but-set-variable \
#    -Wno-unused-parameter"])
AC_PROG_CC
#AC_ADD_CXXFLAGS([-std=gnu++0x]) 
# Checks for libraries.
# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [main])
#ACX_PTHREAD
#LIBS="$PTHREAD_LIBS $LIBS"
#CFLAGS="$CFLAGS $PTHREAD_CFLAGS"
CXXFLAGS="$CXXFLAGS -pthread"
LDFLAGS="$LDFLAGS -pthread"
# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE

# Checks for library functions.
AC_CHECK_FUNCS([pow])

AC_TYPE_SIZE_T

AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_OUTPUT
