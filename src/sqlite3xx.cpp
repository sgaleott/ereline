#include "sqlite3xx.hpp"
#include "datatypes.hpp"
#include "logging.hpp"

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

////////////////////////////////////////////////////////////////////////////////

static int
signalSqlite3ErrorAndThrow(const char * file_name, int code)
{
    auto message = boost::format("Error operating on SQLite3 "
				 "database %1%: \"%2%\"")
	% file_name
	% sqlite3_errstr(code);

    Logger * log = Logger::get_instance();
    log->error(message);

    throw IoError(message);
}

////////////////////////////////////////////////////////////////////////////////

SQLite3Connection::SQLite3Connection(const char * a_file_name)
    : file_name(a_file_name)
{
    Logger * log = Logger::get_instance();

    log->debug(boost::format("Opening connection to SQLite3 database %1%")
	       % a_file_name);
    int result = sqlite3_open(a_file_name, &conn);
    if(result != SQLITE_OK) {
	signalSqlite3ErrorAndThrow(a_file_name, result);
    }

    sqlite3_extended_result_codes(conn, 1);
}

////////////////////////////////////////////////////////////////////////////////

SQLite3Connection::~SQLite3Connection()
{
    sqlite3_close(conn);
}

////////////////////////////////////////////////////////////////////////////////

SQLite3Statement::SQLite3Statement(SQLite3Connection & a_db,
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
	signalSqlite3ErrorAndThrow(db.file_name.c_str(), result);
    }
}

////////////////////////////////////////////////////////////////////////////////

SQLite3Statement::~SQLite3Statement() 
{ 
    Logger * log = Logger::get_instance();
    sqlite3_finalize(ptr);
}

////////////////////////////////////////////////////////////////////////////////

int 
SQLite3Statement::step() 
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
	signalSqlite3ErrorAndThrow(db.file_name.c_str(), result);
    }

    return result;
}
