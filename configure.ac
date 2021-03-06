# file: configure.ac
# author: Olivier Mesnard (mesnardo@gwu.edu)
# brief: Input file for Autoconf.

# set-up
AC_PREREQ([2.69])
AC_INIT([PetIBM-utilities], [0.1dev], [mesnardo@gwu.edu])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_SRCDIR([src/field.cpp])
AC_CONFIG_HEADERS([config/config.h])
AM_SILENT_RULES([no])
AM_MAINTAINER_MODE([disable])
AC_LANG([C++])
LT_PREREQ([2.2])
LT_INIT([dlopen])

# check for programs
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_AWK
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_MKDIR_P
AC_PROG_CC_C_O

# check for typedefs, structures, and compiler characteristics
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_INT64_T
AC_CHECK_TYPES([ptrdiff_t])

# check for header files
AC_CHECK_HEADERS([arpa/inet.h fcntl.h float.h limits.h netdb.h stddef.h stdlib.h string.h strings.h sys/param.h sys/socket.h sys/time.h sys/timeb.h unistd.h wchar.h wctype.h])

# check for library functions
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_FUNC_STRERROR_R
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([dup2 fchdir getcwd getpagesize gettimeofday localtime_r memset mkdir munmap pow regcomp rmdir socket sqrt strcasecmp strchr strdup strerror strrchr strstr strtol strtoull])

AC_SUBST([BUILDDIR], [`pwd`])
CPPFLAGS_="-I$prefix/include"
LDFLAGS_="-L$prefix/lib -Wl,-rpath,$prefix/lib"
CPPFLAGS_PREPEND($CPPFLAGS_)
LDFLAGS_PREPEND($LDFLAGS_)

CONFIGURE_PETSC

# list of Makefiles to generate
AC_CONFIG_FILES([Makefile
                 include/Makefile
                 src/Makefile
                 applications/Makefile
                 applications/convert/Makefile
                 applications/interpolation/Makefile
                 applications/vorticity/Makefile])

# output message
AC_OUTPUT
echo
echo "====================="
echo "Configuration summary"
echo "====================="
AC_MSG_NOTICE([PetIBM-utilities appears to have been configured successfully

Prefix: '${prefix}'

CXX: ${CXX}
CXXFLAGS: ${CXXFLAGS}

CPPFLAGS: ${CPPFLAGS}

LDFLAGS: ${LDFLAGS}

LIBS: ${LIBS}

PETSC_DIR: ${PETSC_DIR}
PETSC_ARCH: ${PETSC_ARCH}

Now type 'make @<:@<target>@:>@'
  where the optional @<:@<target>@:>@ is:
    all            - build all binaries
    install        - install everything
])
