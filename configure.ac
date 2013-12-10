# Process this file with autoconf to produce `configure'

AC_INIT([ereline],[0.1],[tomasi@lambrate.inaf.it])
AC_PREREQ([2.59])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([foreign])
AC_PROG_CXX
AC_PROG_CC

AX_CXX_COMPILE_STDCXX_11([],[mandatory])
AX_BOOST_BASE([1.41])
AX_BOOST_UNIT_TEST_FRAMEWORK

AC_CHECK_LIB(cfitsio, ffopen,, AC_MSG_ERROR(Cannot find the CFITSIO library.))

AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test "$DOXYGEN" == ""; then
   AC_MSG_WARN([Doxygen not found -- the Developer's manual will not be generated])
fi

AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])
AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([doc/developer_manual/Doxyfile])])
	       

# create a config.h file (Automake will add -DHAVE_CONFIG_H)
AC_CONFIG_SRCDIR([src/ereline_test.cpp])
AC_CONFIG_HEADER([config.h])

AC_CONFIG_FILES([Makefile
                 src/Makefile
		 test/Makefile
		 doc/developer_manual/Makefile])
AC_OUTPUT
