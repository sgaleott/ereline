# Process this file with Automake to produce `Makefile.in'

AUTOMAKE_OPTIONS = gnu subdir-objects

bin_PROGRAMS = dx11d_pipeline
dx11d_pipeline_SOURCES = \
	dx11d_pipeline.cpp \
	logging.cpp \
	configuration.cpp \
	ahf_info.cpp \
	data_binning.cpp \
	dipole_parameters.cpp \
	mpi_processes.cpp \
	rotmatrix.cpp \
	planck_velocity.cpp \
	fits_object.cpp \
	dipole_fit.cpp \
	da_capo.cpp \
	mademoiselle.cpp \
	deltavv.cpp \
	gain_table.cpp \
	smooth_gains.cpp \
	project_tod.cpp \
	ringset.cpp \
	io.cpp \
	datatypes.cpp \
	sqlite3xx.cpp \
	squeezer.cpp \
	misc.cpp

AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(SQLITE3_CFLAGS)
LDADD = $(SQLITE3_LDFLAGS) $(BOOST_SYSTEM_LIB) $(BOOST_FILESYSTEM_LIB)
