#include <mpi.h>

#include "data_binning.hpp"

#include "data_binning_results.hpp"
#include "configuration.hpp"
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
#include <limits>

extern "C" {
#include "chealpix.h"
}

/* Return "true" if there were at least two pixels that have been hit
 * during the pointing period.
 *
 * Note: three fields "binned_data.quality_flag" and
 * "binned_data.nside" must already have been set before calling this
 * function! */
static bool
project_pid_into_map(const DifferencedData & diff_data,
                     const PointingData & pointings,
                     const std::vector<double> dipole_tod,
                     const std::vector<double> galactic_tod,
                     const Range_t<size_t> & index_range,
                     Binned_data_t & binned_data)
{
    Logger * log = Logger::get_instance();
    const size_t numPixs = 12 * binned_data.nside * binned_data.nside;

    std::vector<int> hits_map (numPixs);
    std::vector<double> data_map (numPixs);
    std::vector<float> model_map (numPixs);

    // bin the samples and calculate the "binned" dipole
    log->debug(boost::format("Running %s with indexes in "
                             "the range [%d, %d] (there are %d samples "
                             "available). NSIDE is %d")
               % __PRETTY_FUNCTION__
               % index_range.start
               % index_range.end
               % diff_data.sky_load.size()
               % binned_data.nside);

    binned_data.min_dipole = std::numeric_limits<double>::max();
    binned_data.max_dipole = std::numeric_limits<double>::min();
    for (size_t sample_idx = index_range.start;
         sample_idx <= index_range.end;
         sample_idx++)
    {
        if ((diff_data.flags.at(sample_idx) & binned_data.quality_flag) == 0)
        {
            const double dipole = dipole_tod.at(sample_idx);

            long pixel_idx = 0;
            ang2pix_nest(binned_data.nside,
                         pointings.theta.at(sample_idx),
                         pointings.phi.at(sample_idx), 
                         &pixel_idx);

            data_map[pixel_idx] += diff_data.sky_load.at(sample_idx);
            model_map[pixel_idx] += dipole + galactic_tod.at(sample_idx);
            hits_map[pixel_idx]++;

            if (dipole > binned_data.max_dipole)
                binned_data.max_dipole = dipole;
            if (dipole < binned_data.min_dipole)
                binned_data.min_dipole = dipole;

        }
    }

    for (size_t pixel_idx = 0; pixel_idx < numPixs; pixel_idx++)
    {
        if (hits_map[pixel_idx] == 0)
            continue;

        binned_data.pix_index.push_back(pixel_idx);
        binned_data.pix_data_sum.push_back (data_map[pixel_idx]);
        binned_data.pix_model_mean.push_back (model_map[pixel_idx] / 
                                              hits_map[pixel_idx]);
        binned_data.pix_num_of_hits.push_back (hits_map[pixel_idx]);
    }

    log->debug(boost::format("The binning of pID %1% produced a map "
                             "with %2% pixels")
               % binned_data.pointing_id
               % binned_data.pix_index.size());

    return binned_data.pix_index.size() >= 2;
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

inline static std::string
sidelobes_tod_file_path(const Configuration & program_conf,
                        const Lfi_radiometer_t & radiometer,
                        short od)
{
     return
        (boost::format("%s/dipole_fit/tods/sidelobes/%s_sidelobes_OD%04d.fits")
         % program_conf.getWithSubst("common.base_output_dir")
         % radiometer.shortName()
         % od)
        .str();
}

////////////////////////////////////////////////////////////////////////////////

inline static std::string
dipole_tod_file_path(const Configuration & program_conf,
                     const Lfi_radiometer_t & radiometer,
                     short od)
{
    return (boost::format("%s/dipole_fit/tods/dipole/%s_dipole_OD%04d.fits")
         % program_conf.getWithSubst("common.base_output_dir")
         % radiometer.shortName()
         % od)
        .str();
}

////////////////////////////////////////////////////////////////////////////////

static void
save_vector_as_tod(const std::string & file_path,
                   short od,
                   const Lfi_radiometer_t & radiometer,
                   const DifferencedData & datadiff_template,
                   const std::vector<double> vector)
{
    DifferencedData datadiff_to_save = datadiff_template;
    datadiff_to_save.sky_load = vector;
    save_tod(ensure_path_exists(file_path), od, radiometer, datadiff_to_save);
}

////////////////////////////////////////////////////////////////////////////////

std::string
binned_data_file_name(const Configuration & program_conf,
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

/* The argument "pid_range" typically contains all the pointings to be
 * processed by the MPI process, not only the PIDs within the given
 * OD. The implementation of "process_one_od" will silently skip all
 * the PIDs outside "od". */
static std::vector<Binned_data_t>
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
        program_conf.get<uint32_t>("bin_data.quality_flag", 6111248);
    const bool debug_flag = program_conf.get<bool>("bin_data.debug", false);

    const std::string pnt_file_path(pointings_file_path(storage_conf));
    const std::string ddf_file_path(datadiff_file_path(storage_conf));

    PointingData pointings;
    DifferencedData datadiff;

    // Load pointing information and differenced data for this OD
    log->debug(boost::format("Going to read pointings and differenced "
                             "data for OD %1%") % od);
    load_pointings(pnt_file_path, pointings);
    load_differenced_data(ddf_file_path, datadiff);

    if(pointings.obt_time.empty()) {
        log->warning(boost::format("No data for OD %1%, skipping it")
                     % pid_range.start->od);
        return std::vector<Binned_data_t> {};
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
        std::string file_path(sidelobes_tod_file_path(program_conf, radiometer, od));
        save_vector_as_tod(file_path, od, radiometer, datadiff, sidelobes);
    }

    log->debug("Computing the amplitude of the dipole convolved with 4\u03c0 beams");
    std::vector<double> convolved_dipole(
        planck_velocity.getConvolvedDipole(datadiff.scet_time,
                                           pointings.theta,
                                           pointings.phi,
                                           pointings.psi));
    log->debug("...done");
    if(debug_flag) {
        std::string file_path(dipole_tod_file_path(program_conf, radiometer, od));
        save_vector_as_tod(file_path, od, radiometer, datadiff, convolved_dipole);
    }

    // Loop over each pointing period that belongs to the current OD
    std::vector<Binned_data_t> binned_pids;
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

        Binned_data_t binned_data(quality_flag, mask.nside, cur_pid->id);
        if(project_pid_into_map(datadiff, 
                                pointings, 
                                convolved_dipole, 
                                sidelobes, 
                                idx_range, 
                                binned_data))
        {
            binned_pids.push_back(binned_data);

            std::string file_path(binned_data_file_name(program_conf,
                                                        radiometer,
                                                        cur_pid->od,
                                                        cur_pid->id));

            save_binned_data(ensure_path_exists(file_path), radiometer,
                             binned_data);
        } else {
            log->warning(boost::format("Not enough pixels covered at NSIDE %1% "
                                       "for pID %2% (OD %3%), radiometer %4%")
                         % mask.nside
                         % cur_pid->id
                         % cur_pid->od
                         % radiometer.shortName());
        }
    }

    return binned_pids;
}

////////////////////////////////////////////////////////////////////////////////

void run_data_binning(Sqlite_connection_t & ucds,
                      const Lfi_radiometer_t & rad,
                      Configuration & program_conf,
                      Configuration & storage_conf,
                      const std::vector<Pointing_t> & list_of_pointings,
                      Data_binning_results_t & result)
{
    Logger * log = Logger::get_instance();
    const int mpi_size = MPI::COMM_WORLD.Get_size();
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();

    log->info("Starting module binData");
    log->increase_indent();

    Lfi_radiometer_t real_radiometer(
        radiometer_to_use(mpi_rank, rad, program_conf, storage_conf));

    // Load all the inputs needed by this module
    load_map(program_conf.getWithSubst("bin_data.mask"), 1, result.mask);

    ringset galactic_pickup(program_conf.getWithSubst("bin_data.galactic_pickup"),
                            program_conf.get<int>("bin_data.total_convolve_order", 9),
                            false);

    Planck_velocity_t planck_velocity(
        storage_conf.getWithSubst("spacecraft_velocity_file"),
        read_dipole_fit_params(program_conf));
    if(program_conf.get<bool>("dipole_fit.use_pencil_beam", false)) {
        planck_velocity.use_pencil_beam();
    } else {
        load_convolution_params(ucds, real_radiometer, planck_velocity);
    }

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

            result.binned_pids.insert(result.binned_pids.end(),
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

    log->debug(boost::format("Looping binData on the ODs \u2208 [%1%, %2%]"
                             "completed, %3% valid fits found")
               % data_range.od_range.start
               % data_range.od_range.end
               % result.binned_pids.size());

    // TODO: here we should save the binned pids...

    log->decrease_indent();
    log->info("Module binData completed.");
}
