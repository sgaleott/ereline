#ifndef SQLITE3XX_HPP
#define SQLITE3XX_HPP

#include <string>
#include <sqlite3.h>

////////////////////////////////////////////////////////////////////////////////

// Wrapper to sqlite3 in order to allow RAII
struct SQLite3Connection {
    sqlite3 * conn;
    std::string file_name;
    SQLite3Connection(const char * a_file_name);
    ~SQLite3Connection() noexcept;
};

////////////////////////////////////////////////////////////////////////////////

// Wrapper to sqlite3_stmt in order to allow RAII
struct SQLite3Statement {
    SQLite3Connection & db;
    std::string sql;
    sqlite3_stmt * ptr;
    int num_of_rows;

    SQLite3Statement(SQLite3Connection & a_db,
		     const char * a_sql,
		     int num_of_bytes = -1);
    ~SQLite3Statement() noexcept;

    int step();

    int column_count() const { return sqlite3_column_count(ptr); }
    int column_int(int ncol) const { return sqlite3_column_int(ptr, ncol); }
    double column_double(int ncol) const { return sqlite3_column_double(ptr, ncol); }
    std::string column_text(int ncol) const {
	return std::string(reinterpret_cast<const char *>(sqlite3_column_text(ptr, ncol)));
    }
};

#endif
