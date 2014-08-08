#include <mpi.h>

#include "dipole_fit.hpp"

#include "configuration.hpp"
#include "dipole_fit_results.hpp"
#include "healpix_map.hpp"
#include "io.hpp"
#include "logging.hpp"
#include "misc.hpp"
#include "mpi_processes.hpp"
#include "ringset.hpp"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_eigen.h>

#include <algorithm>
#include <fstream>
#include <iostream>

extern "C" {
#include "chealpix.h"
}

bool
fit_model_and_data(const Binned_data_t & binned_data,
                   const std::vector<float> & maskMap,
                   Dipole_fit_t & dipole_fit)
{
  // Compute size of the masked vectors
  size_t maskedLen = 0;
  for (size_t idx = 0; idx < binned_data.pix_model_mean.size(); ++idx) {
      if (maskMap[binned_data.pix_index[idx]] != 0)
          ++maskedLen;
  }

  // Build masked array
  std::vector<double> dipole(maskedLen);
  std::vector<double> data(maskedLen);
  size_t maskedIdx = 0;
  for (size_t idx = 0; idx < binned_data.pix_model_mean.size(); ++idx) {
      if (maskMap[binned_data.pix_index[idx]] != 0)
      {
          dipole[maskedIdx] = binned_data.pix_model_mean[idx];
          data[maskedIdx] = 
              binned_data.pix_data_sum[idx] / binned_data.pix_num_of_hits[idx];
          maskedIdx++;
      }
  }

  // Un-weighted linear fit
  double c0, c1, cov00, cov01, cov11, chisq;
  gsl_fit_linear (dipole.data(), 1,
                  data.data(), 1,
                  maskedLen, &c0, &c1, &cov00, &cov01, &cov11, &chisq);

  // Set gain and offset
  if (c1 > 0) {
      dipole_fit.gainv = c1;
      dipole_fit.offset = c0;
      return true;
  } else {
      dipole_fit.gainv = 0.0;
      dipole_fit.offset = 0.0;
      return false;
  }
}

////////////////////////////////////////////////////////////////////////////////

inline static std::string
gain_table_file_path(const Configuration & program_conf,
                     const Lfi_radiometer_t & radiometer)
{
    return (boost::format("%s/dipole_fit/%s_dipole_fit_gains.fits")
        % program_conf.getWithSubst("common.base_output_dir")
        % radiometer.shortName()).str();
}

////////////////////////////////////////////////////////////////////////////////

std::string
dipole_fit_file_name(const Configuration & program_conf,
                     const Lfi_radiometer_t & radiometer,
                     int od,
                     int pid)
{
    return (boost::format("%s/dipole_fit/fits/%s_fit_OD%04d_pid%06d.fits")
        % program_conf.getWithSubst("common.base_output_dir")
        % radiometer.shortName()
        % od
        % pid)
        .str();
}

////////////////////////////////////////////////////////////////////////////////

static void
extract_gains(const std::vector<Dipole_fit_t> & list_of_fits,
              const Range_t<int> & pid_range,
              Gain_table_t & gain_table)
{
    const size_t num_of_fits = pid_range.end - pid_range.start + 1;
    gain_table.pointingIds.resize(num_of_fits);
    gain_table.gain.resize(num_of_fits);
    gain_table.offset.resize(num_of_fits);

    auto cur_fit = list_of_fits.begin();
    for(size_t idx = 0; idx < num_of_fits; ++idx) {
        const int pid = pid_range.start + idx;
        gain_table.pointingIds[idx] = pid;

        if(cur_fit != list_of_fits.end() && 
           cur_fit->binned_data.pointing_id == pid) {
            gain_table.gain[idx] = 1. / cur_fit->gainv;
            gain_table.offset[idx] = cur_fit->offset;

            ++cur_fit;
        } else {
            gain_table.gain[idx] = 0.0;
            gain_table.offset[idx] = 0.0;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void
run_dipole_fit(Sqlite_connection_t & ucds,
               const Lfi_radiometer_t & rad,
               Configuration & program_conf,
               Configuration & storage_conf,
               const std::vector<Pointing_t> & list_of_pointings,
               const Data_binning_results_t & binned_data,
               Dipole_fit_results_t & result)
{
    Logger * log = Logger::get_instance();
    const int mpi_size = MPI::COMM_WORLD.Get_size();
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();

    log->info("Starting module dipoleFit");
    log->increase_indent();

    Lfi_radiometer_t real_radiometer(
        radiometer_to_use(mpi_rank, rad, program_conf, storage_conf));

    // Copy a few fields from the Data_binning_results_t structure
    // into the Dipole_fit_results_t structure (they'll be handy to Da
    // Capo).
    result.mask = binned_data.mask;
    result.pids_per_process = binned_data.pids_per_process;

    // Run the fitting routine on every pID in binned_data
    result.dipole_fits.clear();
    for(auto const & binned_pid : binned_data.binned_pids) {
        Dipole_fit_t dipole_fit(binned_pid);
        if(fit_model_and_data(binned_pid, 
                              binned_data.mask.pixels,
                              dipole_fit))
        {
            result.dipole_fits.push_back(dipole_fit);
        }
    }

    // So far each MPI process worked on its own pIDs. Now we collect
    // all the gains and offsets in one large vector and save it into
    // a FITS file.

    log->debug("Gathering all gains for this MPI process from the "
               "Dipole_fit_t structures...");
    Data_range_t data_range;
    get_local_data_range(mpi_rank, mpi_size, list_of_pointings,
                         result.pids_per_process, data_range);

    std::vector<Pointing_t>::const_iterator first_pid, last_pid;
    get_pid_iterators_for_range(list_of_pointings, data_range.pid_range,
                                first_pid, last_pid);
    if(first_pid == list_of_pointings.end() ||
       last_pid == list_of_pointings.end() ||
       first_pid->id != data_range.pid_range.start ||
       last_pid->id != data_range.pid_range.end)
    {
        log->error("Mismatch in the pointing IDs");
        abort();
    }

    Gain_table_t gain_table;
    extract_gains(result.dipole_fits,
                  Range_t<int> { first_pid->id, last_pid->id },
                  gain_table);

    log->info("Waiting for all the MPI processes to get here...");
    MPI::COMM_WORLD.Barrier();
    log->info("...done, merging the results");

    // Retrieve the gain table from all the other MPI processes and put
    // them together in results.gain_table
    gain_table.mergeResults();
    // Since odd and even MPI processes work on different radiometer arms
    // (M/S), we discard those gains that do not belong to the arm
    // being analyzed by the current MPI process
    gain_table.selectRadiometerGains(mpi_rank % 2, 2,
                                     result.pids_per_process);

    if(mpi_rank == 0 || mpi_rank == 1) {
        const std::string gain_file_path(gain_table_file_path(program_conf,
                                                              real_radiometer));
        log->info(boost::format("Saving dipoleFit gains into %1%")
                  % gain_file_path);
        save_gain_table(ensure_path_exists(gain_file_path),
                        real_radiometer, gain_table);
    }

    log->decrease_indent();
    log->info("Module dipoleFit completed.");
}
