#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <string>
#include <ostream>
#include <boost/format.hpp>

/** \brief Logging class
 *
 * The Logger class is used to provide informative messages to the
 * user as well as useful debug information to the developer.
 *
 * You do not instantiate objects of the Logger type. Instead, there
 * is a global object that can be accessed using
 * Logger::get_instance().
 *
 * Here is an example which shows how to use it. By default the
 * messages are print to the standard error stream.
 *
 * \code
 * #include <logger.hpp>
 *
 * void do_something()
 * {
 *    auto * log = Logger::get_instance();
 *    log->info("Entering do_something()");
 *
 *    // ...
 *
 *    log->info("Here do_something() ends");
 * }
 * \endcode
 */
class Logger {
public:
    enum class Log_level { DEBUG, INFO, WARNING, ERROR };

    /** \brief Return a pointer to the global Logger object.
     *
     * This static method returns a pointer to the global Logger
     * object, the one and only that must be used. Typically, at the
     * beginning of every function that uses logging facilities you
     * call this function to retrieve a pointer to the object, and
     * then you use it. This object must never be freed.
     */
    static Logger * get_instance();

    /** \brief Emit a log message to the default stream
     *
     * This function is the most generic way to produce log messages.
     * The caller must specify the level of the message (i.e., its
     * "severity") as well as the message itself. The message will be
     * printed only if the default log level is lower or equal to \a
     * level, otherwise the function will do nothing.
     *
     * Typically, the user will use other high-level functions
     * provided by Logger, e.g., Logger::debug, Logger::info,
     * Logger::warning, and Logger::error.
     */
    void log(Log_level level, const std::string & string) const;

    /** \brief Increase the level of indentation for all the following
     * messages
     *
     * This function is useful when the code is entering some section
     * where a lot of error messages are going to be produced. The
     * Logger::log function will then use horizontal spaces to
     * visually suggest a tree-like structure of the messages.
     *
     * To decrease the indentation level, use Logger::decrease_indent().
     *
     * Consider the following example:
     * \code
     * void something_deep()
     * {
     *    log->info("Here I am, within something_deep()");
     *    log->info("And now I quit");
     * }
     *
     * void main()
     * {
     *    auto * log = Logger::get_instance();
     *    log->info("Entering some big block of code...");
     *
     *    log->increase_indent();
     *    something_deep();
     *    log->decrease_indent();
     *
     *    log->info("Back to main");
     * }
     * \endcode
     *
     * The output of the program will be:
     * \verbatim
     2013-01-01T00:00:00   INFO: Entering some big block of code...
     2013-01-01T00:00:00   INFO:     Here I am, within something_deep()
     2013-01-01T00:00:00   INFO:     And now I quit
     2013-01-01T00:00:00   INFO: Back to main
     \endverbatim
     */
    void increase_indent() { indent_level++; }
    /// \brief Decrease the indent level. See Logger::increase_indent().
    void decrease_indent() { if(indent_level > 0) indent_level--; }

    /** \brief Set a different output stream
     *
     * This function allows to select a different destination for
     * logging messages than standard error. You must pass the pointer
     * to a \c std::ostream object whose lifetime is equal to the
     * global Logger object (i.e., probably you want it to be kept
     * open till the program ends).
     *
     * This means that the following code is likely to make the
     * program crash:
     *
     * \code
     * void do_something()
     * {
     *    std::ofstream log_stream("./my_log.txt");
     *    auto * log = Logger::get_instance();
     *    log->set_stream(&log_stream);
     * }
     * \endcode
     *
     * The problem with this code is that \c log_stream will be
     * destroyed once \c do_something() ends: if some other place in
     * the code is going to call Logger's functions after the call to
     * \c do_something(), it will refer to a stream that does no
     * longer exist in memory.
     *
     * The safest place to use Logger::set_stream is in \c main.
     *
     * \note The Logger class will never free the object \a
     * new_log_stream. This means that if you create it using \c new,
     * you must call \c delete before quitting, otherwise the object
     * will never be destroyed (and if it is a buffered file, it will
     * not be flushed).
     */
    void set_stream(std::ostream * new_log_stream) { log_stream = new_log_stream; }

    void debug(const std::string & string) const { log(Log_level::DEBUG, string); }
    void debug(const boost::format & string_fmt) const { debug(boost::str(string_fmt)); }

    void info(const std::string & string) const { log(Log_level::INFO, string); }
    void info(const boost::format & string_fmt) const { info(boost::str(string_fmt)); }

    void warning(const std::string & string) const { log(Log_level::WARNING, string); }
    void warning(const boost::format & string_fmt) const { warning(boost::str(string_fmt)); }

    void error(const std::string & string) const { log(Log_level::ERROR, string); }
    void error(const boost::format & string_fmt) const { error(boost::str(string_fmt)); }

private:
    static bool exist_instance_flag;
    static Logger * singleton;

    std::ostream * log_stream;

    int indent_level;
    Log_level log_level;

    Logger() // Private constructor
	{
	    log_stream = &std::cerr;
	    indent_level = 0;
	    log_level = Log_level::INFO;
	}

public:
    /** \brief Select which log message should be produced and which
     * ones should be ignored.
     *
     * Set the minimum level of severity for the messages. After a
     * call to Logger::set_log_level, every message whose level is
     * smaller than \a new_level will be silently thrown out.
     *
     * This facility is useful if you want to reduce the volume of
     * messages produced by your application.
     *
     * The following example will turn off debug messages before
     * running \c do_something() (presmuably because the function
     * outputs a lot of them), and then it will restore the original
     * log level thanks to Logger::get_log_level():
     * \code
     * Log_level old_level = log->get_log_level();
     * log->set_log_level(Log_level::INFO);
     * do_something();
     * log->set_log_level(old_level);
     * \endcode
     */
    void set_log_level(Log_level new_level);

    Log_level get_log_level() const { return log_level; }
};

#endif
