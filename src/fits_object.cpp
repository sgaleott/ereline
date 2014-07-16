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
  if (overwrite == true && name[0] != '!')
    name.insert(0,"!");
  
  if (fits_create_file(&ptr, const_cast<char*>(name.c_str()), &status))
    fits_report_error(stderr, status);
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



/* 
 * Go to specific HDU 
 */
void 
FitsObject::gotoHDU(const int& hduNumber)
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
    fits_report_error(stderr, status);
  
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
    fits_report_error(stderr, status);
  
  return rows;
}


/*
 * Set keyword (string keyword case)
 */
template <> void 
FitsObject::setKey <std::string> (const std::string& keyName, 
				  const std::string& keyValue, 
				  const std::string& comment) 
{
  int status=0;
  
  fits_update_key (ptr, 
		   fitsType<std::string>(),
		   const_cast <char*> (keyName.c_str()),
		   const_cast <char*> (keyValue.c_str()),
		   const_cast <char*> (comment.c_str()),
		   &status);
  if (status != 0)
    fits_report_error(stderr, status);
}

/*EoF*/
