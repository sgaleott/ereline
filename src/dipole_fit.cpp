#include "dipole_fit.hpp"

#include "configuration.hpp"
#include "healpix_map.hpp"
#include "io.hpp"
#include "logging.hpp"
#include "misc.hpp"
#include "mpi_processes.hpp"
#include "ringset.hpp"
#include "squeezer.hpp"

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

dipoleFit::dipoleFit(uint32_t a_qualityFlag, 
		     int a_nSide, 
		     int a_pointingID)
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
		   const std::vector<uint32_t>& flag,		  
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
	       const std::vector<uint32_t> & flag,
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

static LfiRadiometer
radiometer_to_use(int mpi_rank, 
		  const LfiRadiometer & user_rad,
		  Configuration & program_conf,
		  Configuration & storage_conf)
{
    Logger * log = Logger::get_instance();

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
	real_radiometer = user_rad;
    else
	real_radiometer = user_rad.twinRadiometer();

    setup_variables_for_radiometer(real_radiometer, program_conf);
    setup_variables_for_radiometer(real_radiometer, storage_conf);
    log->info(boost::format("Going to process data for radiometer %1%")
	      % real_radiometer.shortName());

    return real_radiometer;
}

////////////////////////////////////////////////////////////////////////////////

static Data_range_t
get_local_data_range(int mpi_rank, 
		     int mpi_size,
		     const std::vector<Pointing_t> & list_of_pointings)
{
    Logger * log = Logger::get_instance();

    auto list_of_ods = build_od_list(list_of_pointings);

    std::vector<Data_range_t> list_of_data_ranges;
    splitOdsIntoMpiProcesses(mpi_size, list_of_ods, list_of_data_ranges);
    log->info(boost::format("The data to analyze will be split in %1% chunks "
			    "(there are %2% MPI jobs running)")
	      % list_of_data_ranges.size()
	      % mpi_size);

    Data_range_t data_range = list_of_data_ranges.at(mpi_rank / 2);
    log->info(boost::format("Data range for dipoleFit: ODs [%1%, %2%] "
			    "(pointings [%3%, %4%]), number of pointings: %5%")
	      % data_range.od_range.start
	      % data_range.od_range.end
	      % data_range.pid_range.start
	      % data_range.pid_range.end
	      % data_range.num_of_pids);

    return data_range;
}

////////////////////////////////////////////////////////////////////////////////

/* Determine if the pointings and the differenced samples in
 * "pointings" and "datadiff" are compatible or not. (Typically, they
 * are not compatible when they refer to different ODs.) */
static void
assert_consistency(const PointingData & pointings, 
		   const DifferencedData & datadiff)
{
    if(pointings.obt_time.size() != datadiff.obt_time.size()) {
	auto msg = boost::format("Mismatch in the number of pointings and "
				 "differenced samples: %1% against %2%")
	    % pointings.obt_time.size() 
	    % datadiff.obt_time.size();
	throw std::runtime_error(msg.str());
    }

    if(pointings.obt_time.front() != datadiff.obt_time.front()) {
	auto msg = boost::format("OBT times do not match between pointings and "
				 "differenced samples: the first time is "
				 "%1% against %2%")
	    % pointings.obt_time.front() 
	    % datadiff.obt_time.front();
	throw std::runtime_error(msg.str());
    }

    if(pointings.obt_time.back() != datadiff.obt_time.back()) {
	auto msg = boost::format("OBT times do not match between pointings and "
				 "differenced samples: the last time is "
				 "%1% against %2%")
	    % pointings.obt_time.back() 
	    % datadiff.obt_time.back();
	throw std::runtime_error(msg.str());
    }
}

////////////////////////////////////////////////////////////////////////////////

inline std::string
pointings_file_path(const Configuration & storage_conf)
{
    return storage_conf.getWithSubst("pointings.base_path") + "/" +
	storage_conf.getWithSubst("pointings.file_name_mask");
}

////////////////////////////////////////////////////////////////////////////////

inline std::string
datadiff_file_path(const Configuration & storage_conf)
{
    return storage_conf.getWithSubst("differenced_data.base_path") + "/" +
	storage_conf.getWithSubst("differenced_data.file_name_mask");
}

////////////////////////////////////////////////////////////////////////////////

static std::vector<dipoleFit>
process_one_od(const Configuration & program_conf,
	       const Configuration & storage_conf,
	       int od,
	       const Healpix::Map_t<float> & mask,
	       const ringset & galactic_pickup,
	       const PlanckVelocity & planck_velocity,
	       Range_t<std::vector<Pointing_t>::const_iterator> pid_range)
{
    Logger * log = Logger::get_instance();
    const uint32_t quality_flag = 
	program_conf.get<uint32_t>("dipole_fit.quality_flag");

    const std::string pnt_file_path(pointings_file_path(storage_conf));
    const std::string ddf_file_path(datadiff_file_path(storage_conf));

    PointingData pointings;
    DifferencedData datadiff;

    // Load pointing information and differenced data for this OD
    log->debug(boost::format("Going to read pointings and differenced "
			     "data for OD %1%") % od);
    decompress_pointings(pnt_file_path, pointings);
    decompress_differenced_data(ddf_file_path, datadiff);

    if(pointings.obt_time.empty()) {
	log->warning(boost::format("No data for OD %1%, skipping it")
		     % pid_range.start->od);
	return std::vector<dipoleFit> {};
    }

    assert_consistency(pointings, datadiff);
    log->debug(boost::format("Pointings and differenced data "
			     "look consistent. There are %d samples, "
			     "going from OBT time %.0f to %.0f")
	       % pointings.obt_time.size()
	       % pointings.obt_time.front()
	       % pointings.obt_time.back());

    // Loop over each pointing period that belongs to the current OD
    std::vector<dipoleFit> fits;
    for(auto cur_pid = pid_range.start; 
	cur_pid != pid_range.end + 1; 
	cur_pid++) 
    {
	if(cur_pid->od != od)
	    continue;

	log->debug(boost::format("Processing pointing with ID %d "
		       "(OBT times go from %.0f to %0.f)")
		   % cur_pid->id
		   % cur_pid->start_time
		   % cur_pid->end_time);

	PointingData pid_pointings(pointings, 
				   cur_pid->start_time, cur_pid->end_time);
	DifferencedData pid_datadiff(datadiff, 
				     cur_pid->start_time, cur_pid->end_time);

	if(pid_pointings.obt_time.empty()) {
	    log->info("No data for this pointing, skipping it");
	    continue;
	}

	std::vector<double> sidelobes(
	    galactic_pickup.getIntensities(pid_pointings.theta,
					   pid_pointings.phi,
					   pid_pointings.psi));

	std::vector<double> convolved_dipole(
	    planck_velocity.getConvolvedDipole(pid_datadiff.scet_time,
					       pid_pointings.theta,
					       pid_pointings.phi,
					       pid_pointings.psi));

	dipoleFit fitter(quality_flag, mask.nside, cur_pid->id);
	if(fitter.fit(pid_datadiff.sky_load,
		      pid_datadiff.flags,
		      pid_pointings.theta,
		      pid_pointings.phi,
		      convolved_dipole,
		      std::vector<size_t> { 0, pid_datadiff.sky_load.size() - 1 },
		      mask.pixels,
		      sidelobes))
	    fits.push_back(fitter);
    }

    return fits;
}

////////////////////////////////////////////////////////////////////////////////

void
run_dipole_fit(const LfiRadiometer & rad,
	       Configuration & program_conf,
	       Configuration & storage_conf,
	       const std::vector<Pointing_t> & list_of_pointings)
{
    Logger * log = Logger::get_instance();
    const int mpi_size = MPI::COMM_WORLD.Get_size();
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();

    log->info("Starting module dipoleFit");
    log->increase_indent();

    LfiRadiometer real_radiometer(
	radiometer_to_use(mpi_rank, rad, program_conf, storage_conf));

    // Load all the inputs needed by this module
    Healpix::Map_t<float> mask;
    load_map(program_conf.getWithSubst("dipole_fit.mask"), 1, mask);

    ringset galactic_pickup(program_conf.getWithSubst("dipole_fit.galactic_pickup"),
			    program_conf.get<int>("dipole_fit.total_convolve_order"),
			    false);

    PlanckVelocity planck_velocity(storage_conf.getWithSubst("spacecraft_velocity_file"));

    auto data_range(get_local_data_range(mpi_rank, mpi_size, list_of_pointings));
    std::vector<Pointing_t>::const_iterator first_pid, last_pid;
    get_pid_iterators_for_range(list_of_pointings, data_range.pid_range, 
				first_pid, last_pid);
    if(first_pid == list_of_pointings.end() ||
       last_pid == list_of_pointings.end() ||
       first_pid->id != data_range.pid_range.start ||
       last_pid->id != data_range.pid_range.end) 
    {
	log->error("Mismatch in the pointing IDs");
	return;
    }

    /* We should now iterate over each pointing period in the range.
     * The problem is that data are saved in larger chunks, each of
     * them spanning one OD (operational day). Therefore we need two
     * nested loops: the first loops over the range of OD, and within
     * each step loads data one OD long; the second loop iterates over
     * the subset of pointing periods falling within the OD. The inner
     * loop is implemented within the function "process_one_od". */
    for(int od = data_range.od_range.start; od <= data_range.od_range.end; ++od) {
	log->info(boost::format("Processing OD %1%") % od);

	setup_od_variable(od, program_conf);
	setup_od_variable(od, storage_conf);

	// Determine the extents of each pointings within this OD
	try {
	    auto fits =
		process_one_od(program_conf, storage_conf, od,
			       mask, galactic_pickup, planck_velocity,
			       Range_t<std::vector<Pointing_t>::const_iterator> { first_pid, last_pid });
	    if(fits.empty()) {
		log->warning(boost::format("No fits between dipole and TODs "
					   "found, skipping OD %1%")
			     % od);
		continue;
	    }
	}
	catch(std::runtime_error & exc) {
	    log->error(boost::format("Error: %1%. Skipping OD %2%")
		       % exc.what() % od);
	    continue;
	}
	catch(std::bad_alloc & exc) {
	    log->error("I'm not able to allocate memory for pointings"
		       "and/or differenced data, so I'm aborting");
	    break;
	}

	log->info(boost::format("Nothing more to do with OD %1%.") % od);
    }

    log->decrease_indent();
    log->info("Module dipoleFit completed.");
}

