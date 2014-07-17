#ifndef MISC_HPP
#define MISC_HPP

#include <string>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

typedef signed char int8;
typedef unsigned char uint8;
typedef signed short int16;
typedef unsigned short uint16;
typedef unsigned int uint;
typedef signed long long int64;
typedef unsigned long long uint64;

// Colatitude of the solar system motion relative to CMB
// (ecliptical coordinates).
//const double SOLSYSDIR_ECL_THETA=1.765248346;

// Longitude of the solar system motion relative to CMB
// (ecliptical coordinates).
//const double SOLSYSDIR_ECL_PHI=2.995840906;

// Speed of the solar system motion relative to CMB in m/s.
//const double SOLSYSSPEED=369000.0;

/*************************LFI VALUES*****************************/

// Colatitude min= / max=
const double SOLSYSDIR_ECL_THETA=1.7656131194951572;

// Longitude min= / max=
const double SOLSYSDIR_ECL_PHI=2.995889600573578;

// Solar system speed +/- 2
const double SOLSYSSPEED=370082.2332;

/****************************************************************/

// Light speed in m/s (CODATA 2006)
const double SPEED_OF_LIGHT=2.99792458e8;

// Average CMB temperature in K (Fixsen, D. J. 1999, Volume 707, Issue 2, pp. 916-920 (2009))
const double TCMB = 2.72548;

// Boltzmann's constant in J/K (CODATA 2006)
const double K_BOLTZMAN=1.3806504e-23;

// Planck's constant in J s (CODATA 2006)
const double H_PLANCK=6.62606896e-34;

// Healpix empty value
const double HEALPIX_UNDEF=-1.6375e30;

// Pi
const double PI=3.141592653589793238462643383279502884197;
const double HALFPI=1.570796326794896619231321691639751442099;
const double TWOPI=6.283185307179586476925286766559005768394;

// Degree to radiant conversion factor
const double DEG_TO_RAD=PI/180.0;

enum PDT {
  PIPE_INT8    =  0,
  PIPE_UINT8   =  1,
  PIPE_INT16   =  2,
  PIPE_UINT16  =  3,
  PIPE_INT32   =  4,
  PIPE_UINT32  =  5,
  PIPE_INT64   =  6,
  PIPE_UINT64  =  7,
  PIPE_FLOAT32 =  8,
  PIPE_FLOAT64 =  9,
  PIPE_BOOL    = 10,
  PIPE_STRING  = 11,
  PIPE_INVALID = -1 
};

// Returns datatype code
template<typename T> inline PDT pipeType();
template<> inline PDT pipeType<int8>   () { return PIPE_INT8;   }
template<> inline PDT pipeType<uint8>  () { return PIPE_UINT8;  }
template<> inline PDT pipeType<int16>  () { return PIPE_INT16;  }
template<> inline PDT pipeType<uint16> () { return PIPE_UINT16; }
template<> inline PDT pipeType<int>    () { return PIPE_INT32;  }
template<> inline PDT pipeType<uint>   () { return PIPE_UINT32; }
template<> inline PDT pipeType<int64>  () { return PIPE_INT64;  }
template<> inline PDT pipeType<uint64> () { return PIPE_UINT64; }
template<> inline PDT pipeType<float>  () { return PIPE_FLOAT32;}
template<> inline PDT pipeType<double> () { return PIPE_FLOAT64;}
template<> inline PDT pipeType<bool>   () { return PIPE_BOOL;   }
template<> inline PDT pipeType<std::string> () { return PIPE_STRING; }

std::string getDate (void);

std::string trim (const std::string &orig);

std::vector<double> angToCart(double theta, double phi);

std::vector<double> cartToAng(std::vector<double> cart);

std::vector<std::string> getDetectorIds(int frequency);
std::string getDetectorId(int diode);
std::vector<std::string> getAllDetectorIds();
std::vector<std::string> getAllDiodes();
int getDetectorIdasInt(std::string detectorId);

double computeMean(std::vector<double> & input);
double computeVariance(std::vector<double> & input);

/*
 * Template function to transform a number
 * of any type to a string.
 */
template<typename T> std::string dataToString (const T &x)
{
  std::ostringstream strstrm;
  strstrm.precision(20);
  strstrm << x;
  return trim(strstrm.str());
}

std::string intToString(int64 x, unsigned int width);

/*
 * Template function to transform a string
 * to a number of any type
 */
template<typename T> T stringToData (const std::string &x)
{
  std::istringstream strstrm(x);
  T value;
  strstrm >> value;
  return value;
}

/*
 * Merge two vectors. The resulting vector is vector_dest.
 */
template <typename T> inline void fillVector (std::vector<T> & vector_dest, const std::vector<T> &arr_orig)
{
  for (uint i=0; i<arr_orig.size(); ++i)
    vector_dest.push_back(arr_orig[i]);
}

/*
 * Search an element in an ordered array
 */
template <class T>
inline int searchArr(const std::vector<T> &array, T value)
{
  if (array.empty())
    return -1;

  if (value > array[array.size()-1])
    return (static_cast<int>(array.size())-1);

  unsigned int res = static_cast<int>(lower_bound(array.begin(), array.end(), value) - array.begin());
  return res;
}

/*
 * Search an element in an ordered array
 */
template <class T>
inline int searchIdx(const std::vector<T> &array, T value)
{
  if (array.empty())
    return -1;

  if (value > array[array.size()-1])
    return static_cast<int>(array.size())-1;

  unsigned int res = static_cast<unsigned int>(lower_bound(array.begin(), array.end(), value) - array.begin());
  if ((array[res] == value)||(res==0))
    return res;
  return res-1;
}

template <typename T>
inline int linearSearch(std::vector<T> array, T value)
{  
  for (unsigned int idx=0; idx<array.size(); idx++)
    {  
      if (array[idx] == value)
	return idx;
    }
  return -1;
}

std::string getExchangeName (const std::string prefix, const std::string od, 
			     const std::string suffix);

void flagConversion (std::vector<int> inputFlag, std::vector<int> & outFlag, 
		     std::vector<int> & commonFlag, bool isHk=false);

void mpiError (int rankMPI, int start, int stop);

void invert2_eig (std::vector<double> & cc);
std::vector<int> sortAndCount(std::vector<int> & pixels);
void twiddle (unsigned int &v);
unsigned int intRandUniform(unsigned int x, unsigned int y, unsigned int z, unsigned int w);
std::string getTFem (int horn);

template<typename T> std::vector<T> sumVectors (std::vector<T> add1, std::vector<T> add2)
{
  std::vector<T> retVec;
  if (add1.size() != add2.size())
    {
      std::cout << "Error: the vectors must have the same size: " 
		<< add1.size() << " " << add2.size() << std::endl;
      return retVec;
    }

  for (size_t idx=0; idx<add1.size(); ++idx)
    retVec.push_back(add1[idx]+add2[idx]);

  return retVec;
}

template<typename T> std::vector<T> mulVectors (std::vector<T> add1, std::vector<T> add2)
{
  std::vector<T> retVec;
  if (add1.size() != add2.size())
    {
      std::cout << "Error: the vectors must have the same size: "
		<< add1.size() << " " << add2.size() << std::endl;
      return retVec;
    }

  for (size_t idx=0; idx<add1.size(); ++idx)
    retVec.push_back(add1[idx]*add2[idx]);

  return retVec;
}

////////////////////////////////////////////////////////////////////////////////

const std::string &
join_paths(const std::string & a,
	   const std::string & b);

////////////////////////////////////////////////////////////////////////////////

/* Save in "first_idx" and "last_idx" the indexes of the first and
 * last element in "vec" that is between "first" and "last". 
 *
 * If the range defined by [first, last] is not within "vec", both
 * "first_idx" and "last_idx" are initialized to -1.
 *
 * Vector "vec" must be in sorted order. */

template<typename T> void
find_range_of_indexes(std::vector<T> vec, 
		      T first, 
		      T last,
		      int & first_idx, 
		      int & last_idx)
{
    auto start_ptr = std::lower_bound(vec.begin(), vec.end(), first);
    auto end_ptr = std::upper_bound(vec.begin(), vec.end(), last);
    if(start_ptr == vec.end() || end_ptr == vec.end()) {
	first_idx = last_idx = -1;
    } else {
	first_idx = std::distance(vec.begin(), start_ptr);
	if(end_ptr != vec.begin())
	    last_idx = std::distance(vec.begin(), end_ptr) - 1;
	else
	    last_idx = first_idx;
    }
}

#endif
