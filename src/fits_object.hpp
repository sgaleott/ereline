/*
 * Methods to read and write FITS files
 *
 * Copyright (C) 2011 LFI-DPC
 *
 * Author: Daniele Tavagnacco <tavagnacco@oats.inaf.it>
 *
 */

#ifndef FITS_OBJECT_HPP
#define FITS_OBJECT_HPP

#include <vector>
#include <string>
#include <sstream>
#include <typeinfo>
#include <iomanip>
#include <cstring>

#include "fitsio.h"

#include "misc.hpp"

using namespace std;


struct fitscolumn
{
  string name;
  string unit;
  string unitTypeFitsChar;
  int64 repcount;
};

// Returns datatype code
template<typename T> inline int fitsType();
//template<> inline int fitsType <char*>   () { return TSTRING;     }
template<> inline int fitsType <int8>   () { return TBYTE;     }
template<> inline int fitsType <uint8>  () { return TBYTE;     }
template<> inline int fitsType <int16>  () { return TSHORT;    }
template<> inline int fitsType <uint16> () { return TUSHORT;   }
template<> inline int fitsType <int>    () { return TINT;      }
template<> inline int fitsType <uint>   () { return TUINT;     }
template<> inline int fitsType <int64>  () { return TLONGLONG; }
template<> inline int fitsType <uint64> () { return TLONGLONG; }
template<> inline int fitsType <float>  () { return TFLOAT;    }
template<> inline int fitsType <double> () { return TDOUBLE;   } 
template<> inline int fitsType <bool>   () { return TLOGICAL;  }
template<> inline int fitsType <string> () { return TSTRING;   }

template<typename T> inline string fitsTypeC();
template<> inline string fitsTypeC <int8>   () { return "B"; }
template<> inline string fitsTypeC <uint8>  () { return "B"; }
template<> inline string fitsTypeC <int16>  () { return "I"; }
template<> inline string fitsTypeC <uint16> () { return "U"; }
template<> inline string fitsTypeC <int>    () { return "J"; }
template<> inline string fitsTypeC <uint>   () { return "V"; }
template<> inline string fitsTypeC <int64>  () { return "K"; }
template<> inline string fitsTypeC <uint64> () { return "K"; }
template<> inline string fitsTypeC <float>  () { return "E"; }
template<> inline string fitsTypeC <double> () { return "D"; } 
template<> inline string fitsTypeC <bool>   () { return "L"; }
template<> inline string fitsTypeC <string> () { return "A"; }
  
class FitsObject
{
private:
  fitsfile *ptr;
  
public:
  FitsObject(){ptr=NULL;};

  ~FitsObject(){if (ptr != NULL) close();};

  /* create a new fits file */
  void create (const string& fileName, bool overwrite=false);

 
  /* write/update file checksum keyword*/
  void writeChecksum();

  /* open the fits file directly to the first HDU with a table */
  void openTable (const string& name);
  
  /* Close the file */
  void close();

  /* Go to specific HDU */
  void gotoHDU(const int& hduNumber);


  /* Insert a new binary table */
  void insertBINtable (vector<fitscolumn>& columns)
  { insertTable (columns, "xtension", BINARY_TBL); };
  void insertBINtable (vector<fitscolumn>& columns, const string& name)
  { insertTable (columns, name, BINARY_TBL); };


  /* Insert a new ascii table */
  void insertASCtable (vector<fitscolumn>& columns)
  { insertTable (columns, "xtension", ASCII_TBL); };
  
  /* Insert a new table. Create a dummy image HDU if not present in the file */
  void insertTable(vector<fitscolumn>& columns, const string& extname, const int& type);
 
  /* Add a comment in the current HDU header */
  void setComment(const string& comment);

  /* Get the efficent chunck size of the curren table */
  int64 getChunkSize();


  /*   Get the keyName value 
    void getKeyRaw (const string& keyName, void* keyalue, int& fitsType);
  */

  /* Read the data from column in the current HDU */
  template <typename T> void getColumn (const string &columnName, vector<T>& outData, 
					const int64& firstRow, const int64& lastRow);

  template <typename T> void getColumn (const string &columnName, vector<T>& outData)
  { getColumn (columnName, outData, 0, -1); };

  /* Write data in a column in the current HDU */
  template <typename T> void writeColumn (const string& columnName, vector<T>& data, const int64& offset);
  template <typename T> void writeColumn (const int& columnNumber, vector<T>& data, const int64& offset);

  /* Template function to get keyName of any type */
  template <typename T> void getKey(const string& keyName, T& keyValue);
  
  /* Template function to set keyName of any type. */
  template <typename T> void setKey(const string& keyName, const T& keyValue, const string& comment);
  template <typename T> void setKey(const string& keyName, const T& keyValue)
  { setKey(keyName, keyValue, ""); };


};


/* Read the data from column in the current HDU */
template <class T> void 
FitsObject::getColumn (const string &columnName, vector<T>& outData, const int64& firstElem, const int64& lastElem)
{
  int status=0;
  // get column number from name
  int colNum=0;
  if (fits_get_colnum (ptr, CASESEN, const_cast<char*>(columnName.c_str()), &colNum, &status))
    fits_report_error (stderr, status);


  long int startEelem=1;
  long int startRow=1;
  long int nElements = lastElem;


  // if no last row specified, read the entire column
  if(lastElem == -1)
    if (fits_get_num_rows (ptr, &nElements, &status))
      fits_report_error (stderr, status);
  

  int64 repc=1;
  char * tmp1=NULL;
  if (fits_get_bcolparmsll (ptr, colNum, tmp1, tmp1, tmp1, &repc, 0, 0, 0, 0, &status))
    fits_report_error (stderr, status);
  
  if (repc != 1)
    {
      startRow  = firstElem/repc + 1;
      startEelem = firstElem%repc + 1;
      if(lastElem == -1)
	nElements *= repc;
    }
  
  
  // build array
  T * tdata = new T[nElements - firstElem];
  
  // get data in temporary c-array
  if (fits_read_col (ptr, fitsType<T>(), colNum, startRow, startEelem, nElements-firstElem, 0, tdata, 0, &status))
    fits_report_error (stderr, status);
  
  // convert c-array into vector
  outData.assign (tdata, tdata+(nElements-firstElem));
  
  // free memory
  delete [] tdata;

}
  
  



/* Write data in a column in the current HDU */
template <typename T> void 
FitsObject::writeColumn (const string &columnName, vector<T>& data, const int64& offset=1)
{
  int status = 0;
  // get column number from name
  int colNum=0;
  if (fits_get_colnum (ptr, CASESEN, const_cast<char*>(columnName.c_str()), &colNum, &status))
    fits_report_error (stderr, status);
  writeColumn (colNum, data, offset);
}



template <typename T> void 
FitsObject::writeColumn (const int& colNum, vector<T>& data, const int64& offset=1)
{
  int status = 0;
  int firstElem=1;
  
  T* tdata = new T[data.size()];
  for (size_t i=0; i<data.size(); i++)
    tdata[i]=data[i];
  
  if (fits_write_col (ptr, fitsType<T>(), colNum, offset, firstElem, data.size(), tdata, &status))
    fits_report_error (stderr, status);
  
  delete [] tdata;  
}




/* Read a keyword in current HDU header */
template <typename T> void 
FitsObject::getKey(const string& keyName, T& keyValue) 
{
  int status=0;

  char *temp=new char[FLEN_CARD];
  char *comment=NULL;
  fits_read_keyword (ptr, const_cast<char*>(keyName.c_str()), temp, comment, &status);
  if (status != 0)
    fits_report_error(stderr, status);

  // if the key is a string i need to remove ' 
  stringstream tmp (string(temp+1,strlen(temp)-2));
  // if key is not a string, do not remove any char
  if (temp[0] != '\'')
    {
      tmp.str("");
      tmp << temp;
    }

  // conversion
  tmp >> keyValue;
  delete [] temp;
}

/* Write a keyword in current HDU header 
 * - specialized function for string in .cpp file
 */
template <typename T> void 
FitsObject::setKey (const string& keyName, const T& keyValue, const string& comment) 
{
  cout << "Setting key '"<< left << setw(20) << keyName << "' to " << keyValue <<endl;
  int status=0;
  
  T tmpVal;
  stringstream tmp;
  tmp << keyValue;
  tmp >> tmpVal;
  
  fits_update_key (ptr, 
		   fitsType<T>(),
		   const_cast <char*> (keyName.c_str()),
		   &tmpVal,
		   const_cast <char*> (comment.c_str()),
		   &status);
  
  if (status != 0)
    fits_report_error(stderr, status);
}




#endif
/*EoF*/
