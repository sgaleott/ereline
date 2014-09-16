.. Describe how the configuration file is processed

Configuration file
==================


C++ implementation
------------------

The contents of the configuration file in the ``dx11d_pipeline``
program are stored in the :cpp:class:`Configuration` class. The class
used the RapidJSON library (bundled with Ereline) to parse the JSON
file.

.. note::

   An unfortunate design choice of this class is the fact that the
   parameters are extracted from the JSON file "on-demand". This means
   that if a parameter is lacking in the JSON file provided by the
   user, the error will be issued only when the program actually tries
   to access this parameter. A more sensible implementation would
   check that all the necessary parameters are present in the JSON
   file immediately after the program starts, so that errors can be
   reported to the user as soon as possible.



Predefined variables
--------------------

String values in the JSON file can reference other keys in the same
file. There are a number of predefined key/value pairs that are
automatically defined by the program. These are:

================= ================================================= ========
Name              Value                                             Example
================= ================================================= ========
``od``            Number of the current OD                          91
``odNNNN``        As ``od``, but zero-padded to 4 characters        0091
``odNNNNNN``      As ``od``, but zero-padded to 6 characters        000091
``odNNNNNNNN``    As ``od``, but zero-padded to 8 characters        00000091
``horn``          Number of the horn being processed                18
``arm``           One-letter identifier of the current arm          ``M``
``frequency_GHz`` Frequency of the current channel, in GHz          44
``mpi_rank``      Zero-based rank of the current MPI process        2
``mpi_rankNN``    As ``mpi_rank``, but zero-padded to 2 characters  02
``mpi_rankNNNN``  As ``mpi_rank``, but zero-padded to 4 characters  0002
``mpi_size``      Number of MPI processes used in this run          8
``mpi_sizeNN``    As ``mpi_size``, but zero-padded to 2 characters  08
``mpi_sizeNNNN``  As ``mpi_size``, but zero-padded to 4 characters  0008
================= ================================================= ========


Fallbacks
---------

The :cpp:class:`Configuration` class implements "fallbacks": when some
parameters are not found, the class can resort to alternatives. Thus,
for instance, if the key specifying the input data required by some
module is missing, then the class can automatically look for the key
specifying where the output data of the previous module in the
pipeline have been saved. The following callbacks are defined:

============================= =============================
If this key is not defined…   …then this key will be used
============================= =============================
``da_capo.input_gains``       ``dipole_fit.output_gains``
``smooth_gains.input_gains``  ``da_capo.output_gains``
============================= =============================

Internally, a fallback is defined by means of the associative array
``fallbacks``, which is a member of :cpp:class:`Configuration` of type
``std::map<std::string, std::string>``. New fallbacks can be freely
defined after the constructor of :cpp:class:`Configuration` has been
called.
