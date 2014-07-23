#include "sqlite3xx.hpp"
#include "datatypes.hpp"
#include "logging.hpp"

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

////////////////////////////////////////////////////////////////////////////////

static int
signalSqlite3ErrorAndThrow(const char * file_name, sqlite3 * connection)
{
    auto message = boost::format("Error operating on SQLite3 "
                                 "database %1%: \"%2%\"")
        % file_name
        % sqlite3_errmsg(connection);

    Logger * log = Logger::get_instance();
    log->error(message);

    throw IoError(message);
}

////////////////////////////////////////////////////////////////////////////////

Sqlite_connection_t::Sqlite_connection_t(const char * a_file_name)
    : file_name(a_file_name)
{
    Logger * log = Logger::get_instance();

    log->debug(boost::format("Opening connection to SQLite3 database %1%")
               % a_file_name);
    int result = sqlite3_open_v2(a_file_name, &conn, SQLITE_OPEN_READONLY, NULL);
    if(result != SQLITE_OK) {
        signalSqlite3ErrorAndThrow(a_file_name, conn);
    }

    sqlite3_extended_result_codes(conn, 1);
}

////////////////////////////////////////////////////////////////////////////////

Sqlite_connection_t::~Sqlite_connection_t()
{
    sqlite3_close(conn);
}

////////////////////////////////////////////////////////////////////////////////

Sqlite_statement_t::Sqlite_statement_t(Sqlite_connection_t & a_db,
                                       const char * a_sql,
                                       int num_of_bytes)
    : db(a_db), sql(a_sql), num_of_rows(0)
{
    Logger * log = Logger::get_instance();
    log->debug(boost::format("Preparing SQL command \"%1%\" on "
                             "SQLite3 database %2%")
               % sql
               % db.file_name);

    int result = sqlite3_prepare(db.conn,
                                 a_sql,
                                 num_of_bytes,
                                 &ptr,
                                 NULL);
    if(result != SQLITE_OK) {
        signalSqlite3ErrorAndThrow(db.file_name.c_str(), db.conn);
    }
}

////////////////////////////////////////////////////////////////////////////////

Sqlite_statement_t::~Sqlite_statement_t()
{
    Logger * log = Logger::get_instance();
    sqlite3_finalize(ptr);
}

////////////////////////////////////////////////////////////////////////////////

int
Sqlite_statement_t::step()
{
    Logger * log = Logger::get_instance();
    int result = sqlite3_step(ptr);
    switch(result) {
    case SQLITE_ROW: ++num_of_rows; break;
    case SQLITE_DONE: {
        log->debug(boost::format("Command \"%1%\" on SQLite3 "
                                 "database %2% completed, "
                                 "%3% elements were returned")
                   % sql
                   % db.file_name
                   % num_of_rows);
        break;
    }
    default:
        signalSqlite3ErrorAndThrow(db.file_name.c_str(), db.conn);
    }

    return result;
}
