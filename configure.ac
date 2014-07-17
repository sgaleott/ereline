# Process this file with autoconf to produce `configure'

AC_INIT([ereline],[0.1],[tomasi@lambrate.inaf.it])
AC_PREREQ([2.59])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([foreign])
AC_PROG_CXX
AC_PROG_CC

AX_CXX_COMPILE_STDCXX_11
AX_BOOST_BASE([1.41])
if test "$HAVE_BOOST" == ""; then
    AC_MSG_RESULT([BOOST Library found, flags are $BOOST_LDFLAGS])
fi

AX_BOOST_SYSTEM
AX_BOOST_FILESYSTEM
AX_BOOST_UNIT_TEST_FRAMEWORK

AX_LIB_SQLITE3([3.6.0])

AX_PROG_CXX_MPI(true, , [AC_MSG_FAILURE([MPI compiler needed])])

AC_CHECK_LIB(cfitsio, ffopen,, AC_MSG_ERROR(Cannot find the CFITSIO library.))

AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AC_CHECK_LIB([chealpix], [pix2ang_nest])

AC_CHECK_PROGS([DOXYGEN], [doxygen])
if test "$DOXYGEN" == ""; then
   AC_MSG_WARN([Doxygen not found -- the Developer's manual will not be generated])
fi

AM_CONDITIONAL([HAVE_DOXYGEN], [test -n "$DOXYGEN"])
AM_COND_IF([HAVE_DOXYGEN], [AC_CONFIG_FILES([doc/developer_manual/Doxyfile])])
	       

# create a config.h file (Automake will add -DHAVE_CONFIG_H)
AC_CONFIG_SRCDIR([src/dx11d_pipeline.cpp])
AC_CONFIG_HEADER([config.h])

AC_CONFIG_FILES([Makefile
                 src/Makefile
		 test/Makefile
		 doc/developer_manual/Makefile])
AC_OUTPUT
