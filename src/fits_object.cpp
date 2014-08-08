/*
 * Methods to read and write FITS files
 *
 * Copyright (C) 2011 LFI-DPC
 *
 * Author: Daniele Tavagnacco <tavagnacco@oats.inaf.it>
 *
 */

#include "fits_object.hpp"

 /*
  * Create a new fits file
  */
void
FitsObject::create (const std::string& fileName, bool overwrite)
{
  int status=0;
  std::string name=fileName;
  if (overwrite && name[0] != '!')
    name.insert(0,"!");

  if (fits_create_file(&ptr, const_cast<char*>(name.c_str()), &status))
    throw_fits_exception(status);
}

/*
 * Write/update file checksum keyword
 */
void
FitsObject::writeChecksum ()
{
  int status=0;
  if (fits_write_chksum (ptr, &status))
    fits_report_error (stderr, status);
}



/*
 * Open fits file at the first HDU with a table
 */
void
FitsObject::openTable (const std::string& name)
{
  int status=0;
  if (fits_open_table (&ptr, name.c_str(), READONLY, &status))
    fits_report_error (stderr, status);
}




/*
 * Close the current file
 */
void
FitsObject::close()
{
  int status=0;
  if (fits_close_file (ptr, &status))
    fits_report_error (stderr, status);
  ptr=NULL;
}

void
FitsObject::getKey(const std::string& keyName, int & keyValue)
{
    int status = 0;
    if(fits_read_key(ptr, TINT, keyName.c_str(), &keyValue, NULL, &status))
        throw_fits_exception(status);
}

void
FitsObject::getKey(const std::string& keyName, double & keyValue)
{
    int status = 0;
    if(fits_read_key(ptr, TDOUBLE, keyName.c_str(), &keyValue, NULL, &status))
        throw_fits_exception(status);
}

void
FitsObject::getKey(const std::string& keyName, std::string & keyValue)
{
    int status = 0;
    char keyValueAsciiz[FLEN_VALUE + 1];
    if(fits_read_key(ptr, TDOUBLE, keyName.c_str(), &keyValueAsciiz[0],
                     NULL, &status))
        throw_fits_exception(status);

    keyValue = keyValueAsciiz;
}

void
FitsObject::setKey(const std::string& keyName,
                   int keyValue,
                   const std::string& comment)
{
    int status = 0;
    if(fits_update_key(ptr, TINT, keyName.c_str(),
                       &keyValue, comment.c_str(), &status))
        throw_fits_exception(status);
}

void
FitsObject::setKey(const std::string& keyName,
                   double keyValue,
                   const std::string& comment)
{
    int status = 0;
    if(fits_update_key(ptr, TDOUBLE, keyName.c_str(),
                       &keyValue, comment.c_str(), &status))
        throw_fits_exception(status);
}

void
FitsObject::setKey(const std::string& keyName,
                   const std::string & keyValue,
                   const std::string& comment)
{
    int status = 0;
    const char * keyValueAsciiz = keyValue.c_str();
    if(fits_update_key(ptr, TSTRING, keyName.c_str(),
                       const_cast<char*>(keyValueAsciiz),
                       comment.c_str(), &status))
        throw_fits_exception(status);
}

/*
 * Go to specific HDU
 */
void
FitsObject::gotoHDU(int hduNumber)
{
  int status=0;
  int hduType=0;
  if (fits_movabs_hdu (ptr, hduNumber, &hduType, &status))
    fits_report_error (stderr, status);
}



/*
 * Insert a new table
 */
void
FitsObject::insertTable(const std::vector<fitscolumn>& columns,
                        const std::string& extname,
                        int type)
{
  int status=0;

  // string to c-arrays
  char **ttype = new char*[columns.size()];
  char **tunit = new char*[columns.size()];
  char **tform = new char*[columns.size()];

  for (size_t i=0; i<columns.size(); i++)
    {
      ttype[i] = const_cast<char*> (columns[i].name.c_str());
      tunit[i] = const_cast<char*> (columns[i].unit.c_str());
      tform[i] = const_cast<char*> (columns[i].unitTypeFitsChar.c_str());
    }

  // write table
  if (fits_create_tbl (ptr, type, 0, static_cast<int>(columns.size()), ttype, tform, tunit, extname.c_str(), &status))
    fits_report_error (stderr, status);


  // set creation date
  if(fits_write_date (ptr, &status))
    fits_report_error (stderr, status);

  // fix the fitsfile
  writeChecksum();

  // cleanup
  delete [] ttype;
  delete [] tform;
  delete [] tunit;
}


/*
 * Add a comment in the current HDU header
 */
void
FitsObject::setComment(const std::string& comment)
{
  int status=0;

  if (fits_write_comment (ptr, const_cast<char*>(comment.c_str()), &status))
    throw_fits_exception(status);

  // fix the fitsfile
  writeChecksum();
}


/*
 * Get the efficient chunk size to write data in a table
 */
int64
FitsObject::getChunkSize()
{
  int status=0;
  long int rows=0;
  if (fits_get_rowsize(ptr, &rows, &status))
    throw_fits_exception(status);

  return rows;
}
