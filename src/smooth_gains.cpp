#include <mpi.h>

#include "smooth_gains.hpp"
#include "configuration.hpp"
#include "dipole_fit_results.hpp"
#include "da_capo_results.hpp"
#include "smooth_gains_results.hpp"
#include "deltavv.hpp"
#include "io.hpp"
#include "logging.hpp"
#include "misc.hpp"

////////////////////////////////////////////////////////////////////////////////

struct Smoother_parameters_t {
    int hybrid_window_len;
    double hybrid_percent;

    int window_len;
    Range_t<int> window_len_dipole;

    double slow_var_percentile;

    Range_t<double> dipole_range;
};

////////////////////////////////////////////////////////////////////////////////

inline static std::string
dipole_table_file_path(const Configuration & program_conf,
                       const Lfi_radiometer_t & radiometer)
{
    return (boost::format("%s/smoother/%s_dipole_extrema.fits")
        % program_conf.getWithSubst("common.base_output_dir")
        % radiometer.shortName()).str();
}

////////////////////////////////////////////////////////////////////////////////

inline static std::string
gain_table_file_path(const Configuration & program_conf,
                     const Lfi_radiometer_t & radiometer)
{
    return (boost::format("%s/smoother/%s_gain.fits")
        % program_conf.getWithSubst("common.base_output_dir")
        % radiometer.shortName()).str();
}

////////////////////////////////////////////////////////////////////////////////

inline static std::string
slow_gain_table_file_path(const Configuration & program_conf,
                          const Lfi_radiometer_t & radiometer)
{
    return (boost::format("%s/smoother/%s_gain_slow.fits")
        % program_conf.getWithSubst("common.base_output_dir")
        % radiometer.shortName()).str();
}

////////////////////////////////////////////////////////////////////////////////

void
read_smoother_parameters(const Configuration & program_conf,
                         Smoother_parameters_t & params)
{
    params.hybrid_window_len =
        program_conf.get<int>("smooth_gains.fast_variations.window_length");
    params.hybrid_percent =
        program_conf.get<double>("smooth_gains.fast_variations.percent");

    params.window_len_dipole.start =
        program_conf.get<int>("smooth_gains.smooth_window.near_dipole_min");
    params.window_len_dipole.end =
        program_conf.get<int>("smooth_gains.smooth_window.near_dipole_max");
    params.window_len =
        program_conf.get<int>("smooth_gains.smooth_window.window_length");

    params.slow_var_percentile =
        program_conf.get<double>("smooth_gains.slow_var_percentile");

    params.dipole_range.start =
        program_conf.get<double>("smooth_gains.dipole_range.min_value");
    params.dipole_range.end =
        program_conf.get<double>("smooth_gains.dipole_range.max_value");
}

////////////////////////////////////////////////////////////////////////////////

void
load_subsampled_data(Sqlite_connection_t & ucds,
                     const Configuration & program_conf,
                     const Lfi_radiometer_t & radiometer,
                     std::vector<int> & pid,
                     std::vector<double> & ref_data,
                     std::vector<double> & hk_data)
{
    const Range_t<int> pointing_range {
        program_conf.get<int>("smooth_gains.hybrid_fit_first_pid"),
        program_conf.get<int>("smooth_gains.hybrid_fit_last_pid")
        };
    const char * hk_sensor = closest_fp_sensor_to_radiometer(radiometer);
    load_subsampled_ref_and_hk_data(ucds, radiometer, hk_sensor,
                                    pointing_range, pid, ref_data, hk_data);
}

////////////////////////////////////////////////////////////////////////////////

void
run_smooth_gains(Sqlite_connection_t & ucds,
                 const Configuration & program_conf,
                 const Configuration & storage_conf,
                 const Lfi_radiometer_t & rad,
                 const std::vector<Pointing_t> & list_of_pointings,
                 const Dipole_fit_results_t & fit_results,
                 const Da_capo_results_t & da_capo_results,
                 Smooth_gains_results_t & smooth_results)
{
    Logger * log = Logger::get_instance();
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
    const int mpi_size = MPI::COMM_WORLD.Get_size();

    log->info("Starting module smoothGains");
    log->increase_indent();

    Lfi_radiometer_t real_radiometer;
    if(mpi_rank % 2 == 0)
        real_radiometer = rad;
    else
        real_radiometer = rad.twinRadiometer();
    log->info(boost::format("I am going to process radiometer %1%")
              % real_radiometer.shortName());

    // Retrieve the parameters used for the smoothing functions from the
    // configuration file provided by the user
    Smoother_parameters_t smoother_params;
    read_smoother_parameters(program_conf, smoother_params);

    log->info(boost::format("Now I am going to fill the missing fits in "
                            "the raw calibration table with zeroes "
                            "(so far there are %1% pointing periods "
                            "and %2% raw calibration values)")
              % list_of_pointings.size()
              % fit_results.dipole_fits.size());
    Gain_table_t outRawTable;
    Delta_vv_t smoother;
    Gain_table_t outMinMaxTable;
    auto cur_dipole_fit = fit_results.dipole_fits.begin();
    for (auto const & cur_pointing : list_of_pointings)
    {
        const int pid = cur_pointing.id;
        double gain = 0.0;
        double offset = 0.0;
        double dipole = 0.0;
        double min_dipole = 0.0;
        double max_dipole = 0.0;
        if (cur_dipole_fit != fit_results.dipole_fits.end() &&
            cur_dipole_fit->binned_data.pointing_id == pid)
        {
            gain = cur_dipole_fit->gainv;
            offset = cur_dipole_fit->offset;
            dipole = cur_dipole_fit->binned_data.getDipolePeakToPeak();
            min_dipole = cur_dipole_fit->binned_data.min_dipole;
            max_dipole = cur_dipole_fit->binned_data.max_dipole;
            cur_dipole_fit++;
        }

        outRawTable.append({ pid, gain, offset });
        smoother.append(Delta_vv_state_t { pid, gain, dipole });
        outMinMaxTable.append(Gain_state_t { pid, max_dipole, min_dipole });
    }

    MPI::COMM_WORLD.Barrier();

    // Merge Results and Select Gains
    outRawTable.mergeResults();
    outRawTable.selectRadiometerGains(real_radiometer.radiometer, 2,
                                      fit_results.pids_per_process);
    smoother.mergeResults();
    smoother.selectRadiometerGains(real_radiometer.radiometer, 2,
                                   fit_results.pids_per_process);
    outMinMaxTable.mergeResults();
    outMinMaxTable.selectRadiometerGains(real_radiometer.radiometer, 2,
                                         fit_results.pids_per_process);

    if (mpi_rank < 2) {
        save_gain_table(ensure_path_exists(dipole_table_file_path(program_conf, 
                                                                  real_radiometer)),
                        real_radiometer,
                        outMinMaxTable);
    }

    MPI::COMM_WORLD.Barrier();

    /*************************DVV****************************/

    log->info("Computing the hybrid (improved \u03b4V/V) calibration model");
    log->increase_indent();

    // Load Reference Data
    std::vector<int> sub_pid;
    std::vector<double> ref_data, hk_data;
    load_subsampled_data(ucds, program_conf, real_radiometer,
                         sub_pid, ref_data, hk_data);
    log->debug(boost::format("#ref_data = %1%, ref_data[0] = %2%, ref_Data[-1] = %3%")
               % ref_data.size() % ref_data.front() % ref_data.back());
    log->debug(boost::format("#hk_data = %1%, hk_data[0] = %2%, hk_Data[-1] = %3%")
               % hk_data.size() % hk_data.front() % hk_data.back());
    smoother.hybridFit(ref_data, hk_data);

    const std::vector<int> pid(outRawTable.pointingIds);
    const std::vector<double> interpGain(smoother.eval(pid));
    const double meanGain = computeMean(interpGain);

    if(interpGain.size() > 1) {
        log->debug(boost::format("Hybrid gains computed: #interpGain = %1%, "
                                 "interpGain = [%2%, %3%, ..., %4%]")
                   % interpGain.size()
                   % interpGain.at(0)
                   % interpGain.at(1)
                   % interpGain.back());
    }

    MPI::COMM_WORLD.Barrier();
    log->decrease_indent();
    log->info("Hybrid calibration constants computed");

    /********************SAVE DVV GTABLE*********************/

    log->info("Extracting the \u03b4V/V table");
    log->increase_indent();

    Gain_table_t outHybridTable;
    for (auto const & cur_pointing : list_of_pointings)
    {
        double gain = meanGain;
        const int locPid = cur_pointing.id;
        auto pidIter = lower_bound(pid.begin(), pid.end(), locPid);
        if (pidIter != pid.end()) {
            gain = interpGain.at(pidIter - pid.begin());
        }

        outHybridTable.append(Gain_state_t { locPid, gain, 0.0 });
    }
    MPI::COMM_WORLD.Barrier();

    // Save gain table
    outHybridTable.mergeResults();
    outHybridTable.selectRadiometerGains(real_radiometer.radiometer, 2,
                                         fit_results.pids_per_process);

    MPI::COMM_WORLD.Barrier();
    log->decrease_indent();
    log->info(boost::format("\u03b4V/V table computed, %1% values collected")
              % outHybridTable.gain.size());

    /*************************OSGTV***************************/

    log->info("OSGTV section starts here");
    log->increase_indent();

    // Zeroing Hybrid Table for Fast Variations
    log->debug("Calculating the high-frequency part of the curve\u2026");
    const std::vector<double> zeroedHybrid(
        outHybridTable.zeroing(smoother_params.hybrid_window_len,
                               smoother_params.hybrid_percent,
                               smoother.dipole));
    log->debug(boost::format("\u2026done, %1% values computed: [%1%, \u2026, %2%]")
               % zeroedHybrid.size()
               % zeroedHybrid.front()
               % zeroedHybrid.back());

    // Smooth Raw Gains for Slow variations
    log->debug("Smoothing the raw gains\u2026");
    const std::vector<double> smoothedRaw(
        outRawTable.gainSmoothing(smoother_params.window_len_dipole.start,
                                  smoother_params.window_len_dipole.end,
                                  smoother_params.window_len,
                                  smoother_params.slow_var_percentile,
                                  smoother_params.dipole_range.start,
                                  smoother_params.dipole_range.end,
                                  smoother.dipole));
    
    // OR //
    //const std::vector<double> smoothedRaw(
    //  outRaw.Table.gainSmoothing(smoother_params.window_len_dipole.start,
    //				   smoother_params.window_len_dipole.end,
    //				   smoother_params.dipole_range.start,
    //				   smoother_params.dipole_range.end,
    //				   smoother_params.jump_positions
    //				   smoother.dipole));
    // END OR //
    
    log->debug(boost::format("\u2026done, %1% values computed: [%1%, \u2026, %2%]")
               % smoothedRaw.size()
               % smoothedRaw.front()
               % smoothedRaw.back());

    // Smoothed Gains
    log->debug("Summing the high-frequency and low-frequency components\u2026");
    const std::vector<double> smoothedGains =
        sumVectors(smoothedRaw, zeroedHybrid);
    log->debug("\u2026done");

    // Smooth offset
    log->debug("Smoothing the offsets\u2026");
    const std::vector<double> smoothedOffsets(
        outRawTable.offsetSmoothing(smoother_params.window_len_dipole.start,
                                    smoother_params.window_len_dipole.end,
                                    smoother_params.dipole_range.start,
                                    smoother_params.dipole_range.end,
                                    smoother.dipole));
    log->debug(boost::format("\u2026done, %1% values computed: [%1%, \u2026, %2%]")
               % smoothedOffsets.size()
               % smoothedOffsets.front()
               % smoothedOffsets.back());

    MPI::COMM_WORLD.Barrier();
    log->decrease_indent();
    log->info("OSGTV section ends here");

    /****************SAVE GAINS AND REDUCED******************/

    for (auto const & cur_pointing : list_of_pointings)
    {
        const int locPid = cur_pointing.id;
        const int locIdx = outRawTable.getPidIndex (locPid);
        const double gain = smoothedGains[locIdx];
        const double offset = smoothedOffsets[locIdx];

        smooth_results.gain_table.append(Gain_state_t { locPid, gain, offset });

        const double gainSlow = smoothedRaw[locIdx];
        smooth_results.slow_table.append(Gain_state_t { locPid, gainSlow, offset });
    }

    MPI::COMM_WORLD.Barrier();

    // Merge Results and Select Gains
    smooth_results.gain_table.mergeResults();
    smooth_results.gain_table.selectRadiometerGains(real_radiometer.radiometer, 2,
                                                    fit_results.pids_per_process);
    smooth_results.slow_table.mergeResults();
    smooth_results.slow_table.selectRadiometerGains(real_radiometer.radiometer, 2,
                                                    fit_results.pids_per_process);

    if (mpi_rank < 2)
    {
        save_gain_table(ensure_path_exists(gain_table_file_path(program_conf, real_radiometer)),
                        real_radiometer,
                        smooth_results.gain_table);
        save_gain_table(ensure_path_exists(slow_gain_table_file_path(program_conf, real_radiometer)),
                        real_radiometer,
                        smooth_results.slow_table);
    }

    log->decrease_indent();
    log->info("Quitting module smoothGains");
}
