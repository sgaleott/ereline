# Process this file with Automake to produce `Makefile.in'

AUTOMAKE_OPTIONS = gnu subdir-objects

bin_PROGRAMS = ereline_test
ereline_test_SOURCES = \
	ereline_test.cpp \
	logging.cpp \
	configuration.cpp \
	rotmatrix.cpp \
	planck_velocity.cpp \
	fits_object.cpp \
	dipole_fit.cpp \
	da_capo.cpp \
	gain_table.cpp \
	smooth_gains.cpp \
	ringset.cpp \
	io.cpp \
	misc.cpp

AM_CXXFLAGS = $(BOOST_CPPFLAGS)
