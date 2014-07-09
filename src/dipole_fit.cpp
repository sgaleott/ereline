#include "dipole_fit.hpp"
#include "misc.hpp"
#include "configuration.hpp"
#include "logging.hpp"
#include "ringset.hpp"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_eigen.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <mpi.h>

extern "C" {
#include "chealpix.h"
}

dipoleFit::dipoleFit(const int a_qualityFlag, const int a_nSide, const int a_pointingID)
{
  gainv = 0.0;
  offset = 0.0;
  maxDipole = 0.0;
  minDipole = 0.0;
  qualityFlag = a_qualityFlag;
  nSide = a_nSide;
  pointingID = a_pointingID;
}

/*
 * bin the data checkin flags using galactic sidelobes
 */
bool
dipoleFit::binData(const vector<double> & data, const vector<int>& flag,		  
		   const vector<double> & theta, const vector<double> & phi, const vector<double> & dipole, 
		   const vector<size_t> & pidRange, const vector<double> & sidelobes)
{
  int numPixs=12*nSide*nSide;
  vector<int> tmpHits (numPixs,0);
  vector<double> tmpData (numPixs,0);
  vector<float> tmpDipole (numPixs,0);
  vector<double> tmpDipoleConstraint (numPixs,0);

  // bin the samples and calculate the "binned" dipole
  for (size_t sampleNum = pidRange[0]; sampleNum <= pidRange[1]; sampleNum++)
    {
      // select "good" samples
      if ((flag[sampleNum]&qualityFlag) == 0) 
	{
	  // get pixel number in map
	  long pixelNum = 0;  
	  ang2pix_nest(nSide, theta[sampleNum], phi[sampleNum], &pixelNum);
	  
	  tmpData[pixelNum] += data[sampleNum];
	  tmpDipole[pixelNum] += static_cast<float>(dipole[sampleNum]+sidelobes[sampleNum]);
	  
	  if (dipole[sampleNum] > maxDipole)
	    maxDipole = dipole[sampleNum];
	  if (dipole[sampleNum] < minDipole)
	    minDipole = dipole[sampleNum];

	  tmpHits[pixelNum] += 1;
	}
    }

  for (int i=0; i<numPixs; i++)
    {
      if (tmpHits[i] != 0)
	{
	  pixSumData.push_back (tmpData[i]);
	  pixSumDipole.push_back (tmpDipole[i]/static_cast<float>(tmpHits[i]));
	  pixSumHits.push_back (tmpHits[i]);
	  pixIndex.push_back (i);
	}
    }

  if (pixIndex.size() < 2)
    return false;

  return true;
  
}

bool
dipoleFit::fitData(const vector<float> & maskMap)
{ 
  // Compute size of the masked vectors
  size_t maskedLen = 0;
  for (size_t idx=0; idx<pixSumDipole.size(); ++idx)
    if (maskMap[pixIndex[idx]] != 0)
      ++maskedLen;

  // Build masked array
  double * dipole = new double[maskedLen];
  double * data = new double[maskedLen];
  size_t maskedIdx = 0;
  for (size_t idx=0; idx<pixSumDipole.size(); ++idx)
    {
      if (maskMap[pixIndex[idx]] != 0) 
	{
	  dipole[maskedIdx] = static_cast<double>(pixSumDipole[idx]);
	  data[maskedIdx] = pixSumData[idx]/pixSumHits[idx];
	  ++maskedIdx;
	}
    }

  // Un-weighted linear fit
  double c0, c1, cov00, cov01, cov11, chisq;
  gsl_fit_linear (dipole, 1, data, 1, maskedLen, &c0, &c1, &cov00, &cov01, &cov11, &chisq);

  // Set gain and offset
  gainv = c1;
  offset = c0;
  
  // Delete tmp arrays
  delete [] data;
  delete [] dipole;

  if (gainv < 0)
    return false;

  return true;
}

/*
 * fit the dipole using galactic sidelobes
 */
bool
dipoleFit::fit(const vector<double> & data, const vector<int> & flag,
	       const vector<double> & theta, const vector<double> & phi, const vector<double> & dipole, 
	       const vector<size_t> & pidRange, const vector<float> & maskMap, const vector<double> & sidelobes)
{
  if (binData(data, flag, theta, phi, dipole, pidRange, sidelobes))
    return fitData(maskMap);

  return false;
}

void
dipoleFit::setGainV(double a_gainv)
{
  gainv=a_gainv;
}
 
void
dipoleFit::setOffset(double a_offset)
{
  offset=a_offset;
}

void
dipoleFit::setPixSumDipole(const vector<float> & inpArr)
{
  pixSumDipole = inpArr;
}

double 
dipoleFit::getGainV() const
{
  return gainv;
}

double 
dipoleFit::getGain() const
{
  return 1./gainv;
}
 
double 
dipoleFit::getOffset() const
{
  return offset;
}

int 
dipoleFit::getPointingID() const
{
  return pointingID;
}
 
int 
dipoleFit::getNSide() const
{
  return nSide;
}

const vector<int> & 
dipoleFit::getPixIndex() const
{
  return pixIndex;
}

const vector<double> & 
dipoleFit::getPixSumData() const
{
  return pixSumData;
}

const vector<int> & 
dipoleFit::getPixSumHits() const
{
  return pixSumHits;
}

const vector<float> & 
dipoleFit::getPixSumDipole() const
{
  return pixSumDipole;
}

double 
dipoleFit::getDipoleVariance() const
{
  return maxDipole-minDipole;
}

double 
dipoleFit::getMaxDipole() const
{
  return maxDipole;
}

double 
dipoleFit::getMinDipole() const
{
  return minDipole;
}

void
dipoleFit::unload()
{
  vector<double>().swap(pixSumData);
  vector<float>().swap(pixSumDipole);
  vector<int>().swap(pixIndex);
  vector<int>().swap(pixSumHits);
  vector<float>().swap(inputMap);
}

