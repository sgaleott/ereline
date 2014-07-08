Ereline - Planck/LFI calibration pipeline for Erebor
====================================================


Introduction
------------

Ereline is a re-implementation of the Planck/LFI calibration pipeline.
It has been designed for testing different approaches in the
calibration and to checking the robustness of LFI's results.

It is written in C++ and accesses data written in FITS format. The
current plan is to simulate the following features of the instrument:

- Differentiation of noisy total-power data
- Effect of non-ideal beams in the observation of the sky
- Presence of a Doppler effect due to the satellite's motion
- Time-dependent gain instabilities


Dependencies
------------

The code uses C++11 features, so you need a reasonably recent C++
compiler (GCC 4.8.x is ok).

The code needs the Boost library (http://www.boost.org) to compile.
Version 1.41 is the oldest one known to be usable with Ereline. The
following modules are used:

- Boost::date_time
- Boost::format
- Boost::property_tree
- Boost::test

Healpix 3.x must be installed and made available to the "configure"
script. The cleanest way to do this is to specify the CPPFLAGS and
LDFLAGS in the command line used to invoke configure:

    ./configure CPPFLAGS=-I$HEALPIX/include LDFLAGS=-L$HEALPIX/lib

A MPI library is required as well. To use it, specify "CXX=mpic++"
(adjust "mpic++" depending on your MPI Library) in the invocation to
"configure", just as above.

SQLite3 and GNU GSL are two more needed dependency. In both cases, the
"configure" script should be smart enough to find them with no user
intervention.

If you want to build the Developer's manual, you need to install
Doxygen as well (http://www.stack.nl/~dimitri/doxygen).


How to compile the code
-----------------------

Download the code to some directory (e.g. /opt/ereline), then run the
following commands from your home directory:

    $ mkdir ereline_build
    $ cd ereline_build
    $ /opt/ereline/configure \
        --with-boost=PATH_TO_BOOST \
        CPPFLAGS=-I$HEALPIX/include \
        LDFLAGS=-L$HEALPIX/lib CXX=mpic++
    $ make
    $ make check

(note that you do not have to run "configure" and "make" in the same
directory as the source code -- this allows to keep the source in a
read-only folder.)

The last command ("make check") is optional: it ensures that the code
works as expected by running a number of test programs and checking
that their results agrees with the expectations.


How to create the documentation
-------------------------------

The User's manual has been written using Sphinx. You will need to
install Python 2.7 and the appropriate "readthedocs" style. The latter
can be installed using pip:

    $ pip install sphinx_rtd_theme

Then, enter the `doc/user_manual` and issue the command

    $ make html

The result will be placed in the directory `_build/html/'.


Development status
------------------

Currently only a few basic modules have been implemented:

- Logging
- Configuration file reading


License
-------

We have still to decide which license to use. In the long run, it will
probably be an open-source license, but at the moment we want to keep
it private until the whole Planck dataset will be released to the
public.

