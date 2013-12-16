.. highlight:: cpp

.. Developer's guide

Developer's guide
=================

Logging
-------

The source code of Ereline provides the implementation of a C++ class
to write messages to the user. This is done via the
:cpp:class:`Logger` class. It provides a number of useful features:

* It is possible to filter which messages to print (e.g., by avoiding
  the production of debug messages) using a so-called "log level". See
  the :cpp:func:`Logger::set_log_level` method.
* In order to enhance the clarity of the logs for highly structured
  codes, indentation can be set via the
  :cpp:func:`Logger::increase_indent` and
  :cpp:func:`Logger::decrease_indent` methods.
* Messages can be divided into four families:

  1. Debugging messages
  2. Informational messages
  3. Warning messages
  4. Error messages

  The library uses log levels for two purposes:

  1. It prefixes each message with a log level tag
  2. It filters out messages below some importance that has been set
     by the user via the :cpp:func:`Logger::set_log_level`.

The :cpp:class::`Logger` class implements the so-called "singleton
pattern". It only supports one instance of the class at the same time,
and therefore it cannot be instantiated directly. Refer to the class'
documentation to see how to use it in your code.

.. cpp:class:: Logger

  The Logger class is used to provide informative messages to the user
  as well as useful debug information to the developer.

  You do not instantiate objects of the Logger type. Instead, there is
  a global object that can be accessed using Logger::get_instance().

  Here is an example which shows how to use it. By default the
  messages are print to the standard error stream ::

    #include <logger.hpp>

    void do_something()
    {
       autolog = Logger::get_instance();
       log->info("Entering do_something()");

       // ...

       log->info("Here do_something() ends");
    }

.. cpp:function:: Logger * Logger::get_instance()

  This static method returns a pointer to the global Logger object,
  the one and only that must be used. Typically, at the beginning of
  every function that uses logging facilities you call this function
  to retrieve a pointer to the object, and then you use it. This
  object must never be freed.

.. cpp:type:: enum Log_level

  Log levels. They are listed in increasing order of importance. Every
  time you ask the :cpp:class:`Logger` class to emit a message (e.g.,
  via the :cpp:func:`Logger::log` method), the user has to specify the
  level associated with the message. The available levels are:
  ``DEBUG``, ``INFO``, ``WARNING``, and ``ERROR``

.. cpp:function:: void Logger::log(Log_level level, const std::string & string) const

  This function is the most generic way to produce log messages. The
  caller must specify the level of the message (i.e., its "severity")
  as well as the message itself. The message will be printed only if
  the default log level is lower or equal to *level*, otherwise the
  function will do nothing.

  Typically, the user will use other high-level functions provided by
  :cpp:class:`Logger`, e.g., :cpp:func:`Logger::debug`,
  :cpp:func:`Logger::info`, :cpp:func:`Logger::warning`, and
  :cpp:func:`Logger::error`.

.. cpp:function:: void Logger::debug(const std::string & string) const
.. cpp:function:: void Logger::debug(const boost::format & string_fmt) const

  Emit a debug message. This is a wrapper around the
  :cpp:func:`Logger::log` method.

.. cpp:function:: void Logger::info(const std::string & string) const
.. cpp:function:: void Logger::info(const boost::format & string_fmt) const

  Emit an informational message. This is a wrapper around the
  :cpp:func:`Logger::log` method.

.. cpp:function:: void Logger::warning(const std::string & string) const
.. cpp:function:: void Logger::warning(const boost::format & string_fmt) const

  Emit a warning message. This is a wrapper around the
  :cpp:func:`Logger::log` method.

.. cpp:function:: void error(const std::string & string) const
.. cpp:function:: void error(const boost::format & string_fmt) const

  Emit a error message. This is a wrapper around the
  :cpp:func:`Logger::log` method. Such error messages are always
  printed, i.e., irrespectively of the log level set by
  :cpp:func:`Logger::set_log_level`.

.. cpp:function:: void Logger::increase_indent()

  This function is useful when the code is entering some section
  where a lot of error messages are going to be produced. The
  :cpp:func:`Logger::log` function will then use horizontal spaces to
  visually suggest a tree-like structure of the messages.
  
  To decrease the indentation level, use :cpp:func:`Logger::decrease_indent`.

  Consider the following example ::

    void something_deep()
    {
       log->info("Here I am, within something_deep()");
       log->info("And now I quit");
    }
  
    void main()
    {
       auto * log = Logger::get_instance();
       log->info("Entering some big block of code...");
  
       log->increase_indent();
       something_deep();
       log->decrease_indent();
  
       log->info("Back to main");
    }

  
  The output of the program will be:

.. code-block:: none

       2013-01-01T00:00:00   INFO: Entering some big block of code...
       2013-01-01T00:00:00   INFO:     Here I am, within something_deep()
       2013-01-01T00:00:00   INFO:     And now I quit
       2013-01-01T00:00:00   INFO: Back to main

.. cpp:function:: void Logger::decrease_indent()

  Decrease the indent level. See :cpp:func:`Logger::increase_indent`.

.. cpp:function:: void Logger::set_stream(std::ostream * new_log_stream)

  Set a different output stream
 
  This function allows to select a different destination for logging
  messages than standard error. You must pass the pointer to a
  ``std::ostream`` object whose lifetime is equal to the global Logger
  object (i.e., probably you want it to be kept open till the program
  ends).
 
  This means that the following code is likely to make the program
  crash ::
 
    void do_something()
    {
       std::ofstream log_stream("./my_log.txt");
       auto * log = Logger::get_instance();
       log->set_stream(&log_stream);
    }
 
  The problem with this code is that *log_stream* will be destroyed
  once ``do_something()`` ends: if some other place in the code is
  going to call Logger's functions after the call to
  ``do_something()``, it will refer to a stream that does no longer
  exist in memory.
 
  The safest place to use :cpp:func:`Logger::set_stream` is in
  ``main``.
 
.. note::

  The Logger class will never free the object *new_log_stream*. This
  means that if you create it using ``new``, you must call ``delete``
  before quitting, otherwise the object will never be destroyed (and
  if it is a buffered file, it will not be flushed).

.. cpp:function:: void Logger::set_log_level(Log_level new_level)

  Select which log message should be produced and which ones should be
  ignored. After a call to this method, every message whose level is
  smaller than *new_level* will be silently thrown out.
 
  This facility is useful if you want to reduce the volume of messages
  produced by your application.
 
  The following example will turn off debug messages before
  running ``do_something()`` (presmuably because the function
  outputs a lot of them), and then it will restore the original
  log level thanks to :cpp:func:`Logger::get_log_level` ::

    Log_level old_level = log->get_log_level();
    log->set_log_level(Log_level::INFO);
    do_something();
    log->set_log_level(old_level);

  See also :cpp:func:`Logger::get_log_level`.

.. cpp:function:: Log_level Logger::get_log_level() const

  Return the log level set by the previous call to
  :cpp:func:`Logger::set_log_level`.

Configuration files
-------------------

