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

Dipole_fit_t::Dipole_fit_t(uint32_t a_qualityFlag, 
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
Dipole_fit_t::binData(const std::vector<double> & data, 
		      const std::vector<uint32_t>& flag,		  
		      const std::vector<double> & theta, 
		      const std::vector<double> & phi, 
		      const std::vector<double> & dipole, 
		      const Range_t<size_t> & index_range, 
		      const std::vector<double> & sidelobes)
{
  Logger * log = Logger::get_instance();
  const int numPixs=12*nSide*nSide;
  std::vector<int> tmpHits (numPixs,0);
  std::vector<double> tmpData (numPixs,0);
  std::vector<float> tmpDipole (numPixs,0);
  std::vector<double> tmpDipoleConstraint (numPixs,0);

  // bin the samples and calculate the "binned" dipole
  log->debug(boost::format("Running Dipole_fit_t::binData with indexes in "
			   "the range [%d, %d] (there are %d samples "
			   "available). NSIDE is %d")
	     % index_range.start 
	     % index_range.end
	     % data.size()
	     % nSide);
  for (size_t sampleNum = index_range.start; sampleNum <= index_range.end; sampleNum++)
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
Dipole_fit_t::fitData(const std::vector<float> & maskMap)
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
Dipole_fit_t::fit(const std::vector<double> & data, 
		  const std::vector<uint32_t> & flag,
		  const std::vector<double> & theta, 
		  const std::vector<double> & phi, 
		  const std::vector<double> & dipole, 
		  const Range_t<size_t> & index_range, 
		  const std::vector<float> & maskMap, 
		  const std::vector<double> & sidelobes)
{
  if (binData(data, flag, theta, phi, dipole, index_range, sidelobes))
    return fitData(maskMap);

  return false;
}

void
Dipole_fit_t::setPixSumDipole(const std::vector<float> & inpArr)
{
  pixSumDipole = inpArr;
}

double 
Dipole_fit_t::getDipoleVariance() const
{
  return maxDipole-minDipole;
}

void
Dipole_fit_t::unload()
{
  std::vector<double>().swap(pixSumData);
  std::vector<float>().swap(pixSumDipole);
  std::vector<int>().swap(pixIndex);
  std::vector<int>().swap(pixSumHits);
}

////////////////////////////////////////////////////////////////////////////////

static Lfi_radiometer_t
radiometer_to_use(int mpi_rank, 
		  const Lfi_radiometer_t & user_rad,
		  Configuration & program_conf,
		  Configuration & storage_conf)
{
    Logger * log = Logger::get_instance();

    /* The way MPI processes are split in Dipole_fit_t is the following: 
     *
     * 1. Processes with even rank analyze the "main" radiometer;
     *
     * 2. Processes with odd rank analyze the "side" radiometer.
     *
     * (Of course, if the user specified a "side" radiometer in the
     * JSON file, things are reversed.) */
    Lfi_radiometer_t real_radiometer;
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
    log->info(boost::format("Data range for Dipole_fit_t: ODs [%1%, %2%] "
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

/* The argument "pid_range" typically contains all the pointings to be
 * processed by the MPI process, not only the PIDs within the given
 * OD. The implementation of "process_one_od" will silently skip all
 * the PIDs outside "od". */
static std::vector<Dipole_fit_t>
process_one_od(const Configuration & program_conf,
	       const Configuration & storage_conf,
	       const Lfi_radiometer_t & radiometer,
	       int od,
	       const Healpix::Map_t<float> & mask,
	       const ringset & galactic_pickup,
	       const Planck_velocity_t & planck_velocity,
	       Range_t<std::vector<Pointing_t>::const_iterator> pid_range)
{
    Logger * log = Logger::get_instance();
    const uint32_t quality_flag = 
	program_conf.get<uint32_t>("dipole_fit.quality_flag");
    const bool debug_flag = program_conf.get<bool>("dipole_fit.debug");

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
	return std::vector<Dipole_fit_t> {};
    }

    assert_consistency(pointings, datadiff);
    log->debug(boost::format("Pointings and differenced data "
			     "look consistent. There are %d samples, "
			     "going from OBT time %.0f to %.0f")
	       % pointings.obt_time.size()
	       % pointings.obt_time.front()
	       % pointings.obt_time.back());

    log->debug("Computing Galactic pickup through sidelobes...");
    std::vector<double> sidelobes(
	galactic_pickup.getIntensities(pointings.theta,
				       pointings.phi,
				       pointings.psi));
    log->debug("...done");
    if(debug_flag) {
	std::string file_path = 
	    (boost::format("%s/dipole_fit/tods/sidelobes/%s_sidelobes_OD%04d.fits")
	     % program_conf.getWithSubst("common.base_output_dir")
	     % radiometer.shortName()
	     % od)
	    .str();
	DifferencedData galactic_datadiff;
	galactic_datadiff.obt_time = datadiff.obt_time;
	galactic_datadiff.scet_time = datadiff.scet_time;
	galactic_datadiff.sky_load = sidelobes;
	galactic_datadiff.flags = datadiff.flags;
	save_tod(ensure_path_exists(file_path), od, radiometer, galactic_datadiff);
    }

    log->debug("Computing the amplitude of the dipole convolved with 4\u03c0 beams");
    std::vector<double> convolved_dipole(
	planck_velocity.getConvolvedDipole(datadiff.scet_time,
					   pointings.theta,
					   pointings.phi,
					   pointings.psi));
    log->debug("...done");
    if(debug_flag) {
	std::string file_path = 
	    (boost::format("%s/dipole_fit/tods/dipole/%s_dipole_OD%04d.fits")
	     % program_conf.getWithSubst("common.base_output_dir")
	     % radiometer.shortName()
	     % od)
	    .str();
	DifferencedData dipole_datadiff;
	dipole_datadiff.obt_time = datadiff.obt_time;
	dipole_datadiff.scet_time = datadiff.scet_time;
	dipole_datadiff.sky_load = convolved_dipole;
	dipole_datadiff.flags = datadiff.flags;
	save_tod(ensure_path_exists(file_path), od, radiometer, dipole_datadiff);
    }

    // Loop over each pointing period that belongs to the current OD
    std::vector<Dipole_fit_t> fits;
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

	const Range_t<uint64_t> obt_range { 
	    cur_pid->start_time, cur_pid->end_time };
	auto idx_range =
	    find_boundaries_in_obt_times(pointings.obt_time, obt_range);

	Dipole_fit_t fitter(quality_flag, mask.nside, cur_pid->id);
	if(fitter.fit(datadiff.sky_load,
		      datadiff.flags,
		      pointings.theta,
		      pointings.phi,
		      convolved_dipole,
		      idx_range,
		      mask.pixels,
		      sidelobes)) {
	    fits.push_back(fitter);
	    if(debug_flag) {
		std::string file_path = 
		    (boost::format("%s/dipole_fit/fits/%s_fit_OD%04d_pid%06d.fits")
		     % program_conf.getWithSubst("common.base_output_dir")
		     % radiometer.shortName()
		     % cur_pid->od
		     % cur_pid->id)
		    .str();
#
		save_dipole_fit(ensure_path_exists(file_path), radiometer,
				fitter);
	    }
	}
    }

    return fits;
}

////////////////////////////////////////////////////////////////////////////////

static void
extract_gains(const std::vector<Dipole_fit_t> & list_of_fits, 
	      Gain_table_t & gain_table)
{
    const size_t num_of_fits = list_of_fits.size();
    gain_table.pointingIds.resize(num_of_fits);
    gain_table.gain.resize(num_of_fits);
    gain_table.offset.resize(num_of_fits);

    for(size_t idx = 0; idx < num_of_fits; ++idx) {
	const auto & cur_fit = list_of_fits.at(idx);

	gain_table.pointingIds[idx] = cur_fit.pointingID;
	gain_table.gain[idx] = 1. / cur_fit.gainv;
	gain_table.offset[idx] = cur_fit.offset;
    }
}

////////////////////////////////////////////////////////////////////////////////

void
run_dipole_fit(Sqlite_connection_t & ucds,
	       const Lfi_radiometer_t & rad,
	       Configuration & program_conf,
	       Configuration & storage_conf,
	       const std::vector<Pointing_t> & list_of_pointings,
	       Dipole_fit_results_t & result)
{
    Logger * log = Logger::get_instance();
    const int mpi_size = MPI::COMM_WORLD.Get_size();
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();

    if(! program_conf.get<bool>("dipole_fit.run")) {
	log->info("dipoleFit will not be run");
	return;
    }

    log->info("Starting module dipoleFit");
    log->increase_indent();

    Lfi_radiometer_t real_radiometer(
	radiometer_to_use(mpi_rank, rad, program_conf, storage_conf));

    // Load all the inputs needed by this module
    load_map(program_conf.getWithSubst("dipole_fit.mask"), 1, result.mask);

    ringset galactic_pickup(program_conf.getWithSubst("dipole_fit.galactic_pickup"),
			    program_conf.get<int>("dipole_fit.total_convolve_order"),
			    false);

    Planck_velocity_t planck_velocity(
	storage_conf.getWithSubst("spacecraft_velocity_file"),
	read_dipole_fit_params(program_conf));
    load_convolution_params(ucds, real_radiometer, planck_velocity);

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
	    auto od_fits =
		process_one_od(program_conf, storage_conf, real_radiometer, od,
			       result.mask, galactic_pickup, planck_velocity,
			       Range_t<std::vector<Pointing_t>::const_iterator> { first_pid, last_pid });
	    if(od_fits.empty()) {
		log->warning(boost::format("No fits between dipole and TODs "
					   "found, skipping OD %1%")
			     % od);
		continue;
	    }

	    result.list_of_fits.insert(result.list_of_fits.end(), 
				       od_fits.begin(), od_fits.end());
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

    extract_gains(result.list_of_fits, result.gain_table);
    log->decrease_indent();
    log->info("Module dipoleFit completed.");
}

