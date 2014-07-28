#ifndef MISC_HPP
#define MISC_HPP

#include <string>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

#include <boost/format.hpp>

typedef signed char int8;
typedef unsigned char uint8;
typedef signed short int16;
typedef unsigned short uint16;
typedef unsigned int uint;
typedef signed long long int64;
typedef unsigned long long uint64;

/****************************************************************/

// Light speed in m/s (CODATA 2006)
const double SPEED_OF_LIGHT=2.99792458e8;

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

std::string trim (const std::string &orig);

void angToCart(double theta, double phi, double result[3]);

std::vector<double> cartToAng(std::vector<double> cart);

std::vector<std::string> getDetectorIds(int frequency);
std::string getDetectorId(int diode);
int getDetectorIdasInt(const std::string & detectorId);

double computeMean(const std::vector<double> & input);
double computeVariance(const std::vector<double> & input);

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

void mpiError (int rankMPI, int start, int stop);

void invert2_eig (std::vector<double> & cc);
std::vector<int> sortAndCount(std::vector<int> & pixels);

struct Lfi_radiometer_t;
const char * closest_fp_sensor_to_radiometer (const Lfi_radiometer_t & rad);

template<typename T> std::vector<T> sumVectors (std::vector<T> add1, std::vector<T> add2)
{
  if (add1.size() != add2.size())
      throw std::runtime_error(boost::format("Error: the vectors must have the same size: %d / %d")
                               % add1.size() % add2.size());

  std::vector<T> retVec(add1.size());
  for (size_t idx = 0; idx < retVec.size(); ++idx)
    retVec[idx] = add1[idx] + add2[idx];

  return retVec;
}

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
