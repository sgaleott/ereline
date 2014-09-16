.. How to use the logging facilities provided by ``dx11d_pipeline``

Logging facilities
==================

Each Ereline module uses logging functions to monitor the
computations. Each log message is a text string associated with the
date and time when the text was produced, and a "severity" level. The
logging library used by Ereline has the following features:
* The user can select which severity levels should be produced and
which ones should be discarded;
* Log messages can be redirected to a file;
* Each MPI process can have its own log file.

In the ``common`` section of the parameter file, the user can select
the minimum severity level of the log messages that need to be
produced by means of the ``log_level`` key. Its value must be a number
between 1 and 4:

===== ==============================================================
Value What is reported to the user
===== ==============================================================
1     Only errors
2     Errors and warnings
3     Errors, warnings, and general information
4     Errors, warnings, general information, and debugging messages
===== ==============================================================

Logging messages in MPI processes
---------------------------------

The user can select a path to a log file. The variable ``{mpi_rank}``
should be used, either in the path or in the file name, in order to
prevent multiple MPI processes from overwriting each other's file.

C++ usage
---------

The :cpp:class:`Logger` class implements the `Singleton pattern
<http://en.wikipedia.org/wiki/Singleton_pattern>`_, in order to have
one and only one instance of this class available during the execution
of the code. The programmer can access such instance by means of
:cpp:func:`Logger::get_instance`:

.. code-block:: c++

   #include "logging.hpp"

   void test()
   {
       Logger * log = Logger::get_instance();
       log->debug("Now I am into the test() function");
   }

When entering a function which is expected to produce plenty of log
messages, it is advisable to mark such messages with an increased
indent level. This is the purpose of the methods
:cpp:func:`Logger::increase_indent` and
:cpp:func:`Logger::decrease_indent`.

.. cpp:function:: void Logger::increase_indent()

   Increase the indent level of the following logging messages by some
   constant amount (typically, 4 spaces).

.. cpp:function:: void Logger::decrease_indent()

   Decrease the indent level of the following logging messages by some
   constant amount (typically, 4 spaces). If no indent was present,
   this is a no-op.

Here is an example:

.. code-block:: c++

   void very_complex_fn()
   {
       Logger * log = Logger::get_instance();
       log->info("Starting a very long and complex process…");
       log->increase_indent();

       log->debug("Doing step 1…");
       // Lots of stuff here
       // ...

       log->debug("Doing step 2…");
       // Lots of stuff here too
       // ...

       log->debug("Done!");
       log->decrease_indent();
   }
