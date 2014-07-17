#include "dipole_fit.hpp"

#include "configuration.hpp"
#include "healpix_map.hpp"
#include "io.hpp"
#include "logging.hpp"
#include "misc.hpp"
#include "mpi_processes.hpp"
#include "ringset.hpp"
#include "squeezer.hpp"

#include <boost/filesystem.hpp>

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
dipoleFit::binData(const std::vector<double> & data, 
		   const std::vector<int>& flag,		  
		   const std::vector<double> & theta, 
		   const std::vector<double> & phi, 
		   const std::vector<double> & dipole, 
		   const std::vector<size_t> & pidRange, 
		   const std::vector<double> & sidelobes)
{
  const int numPixs=12*nSide*nSide;
  std::vector<int> tmpHits (numPixs,0);
  std::vector<double> tmpData (numPixs,0);
  std::vector<float> tmpDipole (numPixs,0);
  std::vector<double> tmpDipoleConstraint (numPixs,0);

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
dipoleFit::fitData(const std::vector<float> & maskMap)
{ 
  // Compute size of the masked vectors
  size_t maskedLen = 0;
  for (size_t idx=0; idx<pixSumDipole.size(); ++idx)
    if (maskMap[pixIndex[idx]] != 0)
      ++maskedLen;

  // Build masked array
  std::vector<double> dipole(maskedLen);
  std::vector<double> data(maskedLen);
  size_t maskedIdx = 0;
  for (size_t idx=0; idx<pixSumDipole.size(); ++idx)
    {
      if (maskMap[pixIndex[idx]] != 0) 
	{
	  dipole[maskedIdx] = pixSumDipole[idx];
	  data[maskedIdx] = pixSumData[idx] / pixSumHits[idx];
	  ++maskedIdx;
	}
    }

  // Un-weighted linear fit
  double c0, c1, cov00, cov01, cov11, chisq;
  gsl_fit_linear (dipole.data(), 1, 
		  data.data(), 1, 
		  maskedLen, &c0, &c1, &cov00, &cov01, &cov11, &chisq);

  // Set gain and offset
  gainv = c1;
  offset = c0;

  if (gainv < 0)
    return false;

  return true;
}

/*
 * fit the dipole using galactic sidelobes
 */
bool
dipoleFit::fit(const std::vector<double> & data, 
	       const std::vector<int> & flag,
	       const std::vector<double> & theta, 
	       const std::vector<double> & phi, 
	       const std::vector<double> & dipole, 
	       const std::vector<size_t> & pidRange, 
	       const std::vector<float> & maskMap, 
	       const std::vector<double> & sidelobes)
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
dipoleFit::setPixSumDipole(const std::vector<float> & inpArr)
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

const std::vector<int> & 
dipoleFit::getPixIndex() const
{
  return pixIndex;
}

const std::vector<double> & 
dipoleFit::getPixSumData() const
{
  return pixSumData;
}

const std::vector<int> & 
dipoleFit::getPixSumHits() const
{
  return pixSumHits;
}

const std::vector<float> & 
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
  std::vector<double>().swap(pixSumData);
  std::vector<float>().swap(pixSumDipole);
  std::vector<int>().swap(pixIndex);
  std::vector<int>().swap(pixSumHits);
  std::vector<float>().swap(inputMap);
}

////////////////////////////////////////////////////////////////////////////////

/* Important note: variables like "horn" and "arm" need to be already
 * set in "program_conf" before calling "run_dipole_fit". */
void
run_dipole_fit(const LfiRadiometer & rad,
	       Configuration & program_conf,
	       Configuration & storage_conf,
	       const std::vector<Od_t> list_of_ods)
{
    Logger * log = Logger::get_instance();
    const int mpi_size = MPI::COMM_WORLD.Get_size();
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();

    log->info("Starting module dipoleFit");

    /* The way MPI processes are split in dipoleFit is the following: 
     *
     * 1. Processes with even rank analyze the "main" radiometer;
     *
     * 2. Processes with odd rank analyze the "side" radiometer.
     *
     * (Of course, if the user specified a "side" radiometer in the
     * JSON file, things are reversed.) */
    LfiRadiometer real_radiometer;
    if(mpi_rank % 2 == 0)
	real_radiometer = rad;
    else
	real_radiometer = rad.twinRadiometer();
    setup_variables_for_radiometer(real_radiometer, program_conf);
    setup_variables_for_radiometer(real_radiometer, storage_conf);

    Healpix::Map_t<int> mask;
    load_map(program_conf.getWithSubst("dipole_fit.mask"), 1, mask);

    ringset galactic_pickup(program_conf.getWithSubst("dipole_fit.galactic_pickup"),
			    program_conf.get<int>("dipole_fit.total_convolve_order"),
			    false);

    std::vector<Data_range_t> list_of_data_ranges;
    splitOdsIntoMpiProcesses(mpi_size, list_of_ods, list_of_data_ranges);
    log->info(boost::format("The data to analyze will be split in %1% chunks "
			    "(there are %2% MPI jobs running)")
	      % list_of_data_ranges.size()
	      % mpi_size);

    const Data_range_t data_range = list_of_data_ranges.at(mpi_rank);
    log->info(boost::format("Range of data for dipoleFit: ODs [%1%, %2%] "
			    "(pointings [%3%, %4%]), number of pointings: %5%")
	      % data_range.od_range.start
	      % data_range.od_range.end
	      % data_range.pid_range.start
	      % data_range.pid_range.end
	      % data_range.num_of_pids);

    for(int od = data_range.od_range.start; od <= data_range.od_range.end; ++od) {
	setup_od_variable(od, program_conf);
	setup_od_variable(od, storage_conf);

	boost::filesystem::path file_path(storage_conf.getWithSubst("pointings.base_path"));
	file_path /= storage_conf.getWithSubst("pointings.file_name_mask");
	PointingData pointings;
	log->info(boost::format("Loading pointings for OD %1% from file %2%")
		  % od % file_path);
	try {
	    decompress_pointings(file_path.string(), pointings);
	}
	catch(SqueezerError & err) {
	    log->error(err.what());
	    continue;
	}

	//DifferencedData datadiff;
	//decompress_differenced_data();
    }

    log->info("Quitting module dipoleFit");
}

