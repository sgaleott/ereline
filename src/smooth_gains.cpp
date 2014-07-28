#include "smooth_gains.hpp"
#include "configuration.hpp"
#include "dipole_fit_results.hpp"
#include "da_capo_results.hpp"
#include "deltavv.hpp"
#include "io.hpp"
#include "logging.hpp"
#include "misc.hpp"

#include <mpi.h>

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
        program_conf.get<int>("smooth_gains.smooth_window.dipole_min");
    params.window_len_dipole.end =
        program_conf.get<int>("smooth_gains.smooth_window.dipole_max");
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
                 const Da_capo_results_t & da_capo_results)
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

    // Retrieve the parameters used for the smoothing functions from the
    // configuration file provided by the user
    Smoother_parameters_t smoother_params;
    read_smoother_parameters(program_conf, smoother_params);

    Gain_table_t outRawTable;
    Delta_vv_t smoother;
    Gain_table_t outMinMaxTable;
    auto cur_dipole_fit = fit_results.list_of_fits.begin();
    for (auto const & cur_pointing : list_of_pointings)
    {
        const int pid = cur_pointing.id;
        double gain = 0.0;
        double offset = 0.0;
        double dipole = 0.0;
        double minDipole = 0.0;
        double maxDipole = 0.0;
        if (cur_dipole_fit != fit_results.list_of_fits.end() &&
            cur_dipole_fit->pointingID == pid)
        {
            gain = cur_dipole_fit->gainv;
            offset = cur_dipole_fit->offset;
            dipole = cur_dipole_fit->getDipoleVariance();
            minDipole = cur_dipole_fit->minDipole;
            maxDipole = cur_dipole_fit->maxDipole;
            cur_dipole_fit++;
        }

        outRawTable.appendPidValue(pid);
        outRawTable.appendGainValue(gain);
        outRawTable.appendOffsetValue(offset);
        smoother.appendPidValue(pid);
        smoother.appendGainValue(gain);
        smoother.appendDipoleValue(dipole);
        outMinMaxTable.appendPidValue(pid);
        outMinMaxTable.appendGainValue(maxDipole);
        outMinMaxTable.appendOffsetValue(minDipole);
    }

    MPI::COMM_WORLD.Barrier();

    // Merge Results and Select Gains
    outRawTable.mergeResults();
    outRawTable.selectRadiometerGains(mpi_rank % 2, 2, fit_results.pids_per_process);
    smoother.mergeResults();
    smoother.selectRadiometerGains(mpi_rank % 2, 2, fit_results.pids_per_process);
    outMinMaxTable.mergeResults();
    outMinMaxTable.selectRadiometerGains(mpi_rank % 2, 2, fit_results.pids_per_process);

    if (mpi_rank < 2) {
        save_gain_table(ensure_path_exists(dipole_table_file_path(program_conf, real_radiometer)),
                        real_radiometer,
                        outMinMaxTable);
    }

    MPI::COMM_WORLD.Barrier();

    /*************************DVV****************************/

    // Load Reference Data
    std::vector<int> sub_pid;
    std::vector<double> ref_data, hk_data;
    load_subsampled_data(ucds, program_conf, real_radiometer,
                         sub_pid, ref_data, hk_data);
    smoother.hybridFit(ref_data, hk_data);

    std::vector<int> pid = outRawTable.pointingIds;
    std::vector<double> interpGain = smoother.eval(pid);
    const double meanGain = computeMean(interpGain);

    MPI::COMM_WORLD.Barrier();

    /********************SAVE DVV GTABLE*********************/

    Gain_table_t outHybridTable;
    for (auto const & cur_pointing : list_of_pointings)
    {
        double gain = meanGain;
        const int locPid = cur_pointing.id;
        auto pidIter = lower_bound(pid.begin(), pid.end(), locPid);
        if (pidIter != pid.end())
            gain = interpGain.at(pidIter - pid.begin());

        outHybridTable.appendPidValue(locPid);
        outHybridTable.appendGainValue(gain);
        outHybridTable.appendOffsetValue(0.0);
    }

    MPI::COMM_WORLD.Barrier();

    // Save gain table
    outHybridTable.mergeResults();
    outHybridTable.selectRadiometerGains(mpi_rank % 2, 2,
                                         fit_results.pids_per_process);

    MPI::COMM_WORLD.Barrier();

    /*************************OSGTV***************************/

    // Zeroing Hybrid Table for Fast Variations
    const std::vector<double> zeroedHybrid(
        outHybridTable.zeroing(smoother_params.hybrid_window_len,
                               smoother_params.hybrid_percent,
                               smoother.dipole));

    // Smooth Raw Gains for Slow variations
    const std::vector<double> smoothedRaw(
        outRawTable.gainSmoothing(smoother_params.window_len_dipole.start,
                                  smoother_params.window_len_dipole.end,
                                  smoother_params.window_len,
                                  smoother_params.slow_var_percentile,
                                  smoother_params.dipole_range.start,
                                  smoother_params.dipole_range.end,
                                  smoother.dipole));

    // Smoothed Gains
    const std::vector<double> smoothedGains =
        sumVectors(smoothedRaw, zeroedHybrid);

    // Smooth offset
    const std::vector<double> smoothedOffsets(
        outRawTable.offsetSmoothing(smoother_params.window_len_dipole.start,
                                    smoother_params.window_len_dipole.end,
                                    smoother_params.dipole_range.start,
                                    smoother_params.dipole_range.end,
                                    smoother.dipole));

    MPI::COMM_WORLD.Barrier();

    /****************SAVE GAINS AND REDUCED******************/

    Gain_table_t gTable;
    Gain_table_t slowTable;
    for (auto const & cur_pointing : list_of_pointings)
    {
        const int locPid = cur_pointing.id;
        const int locIdx = outRawTable.getPidIndex (locPid);
        const double gain = smoothedGains[locIdx];
        const double offset = smoothedOffsets[locIdx];

        gTable.appendPidValue (locPid);
        gTable.appendGainValue (gain);
        gTable.appendOffsetValue (offset);

        const double gainSlow = smoothedRaw[locIdx];

        slowTable.appendPidValue (locPid);
        slowTable.appendGainValue (gainSlow);
        slowTable.appendOffsetValue (offset);
    }

    MPI::COMM_WORLD.Barrier();

    // Merge Results and Select Gains
    gTable.mergeResults();
    gTable.selectRadiometerGains(mpi_rank % 2, 2, fit_results.pids_per_process);
    slowTable.mergeResults();
    slowTable.selectRadiometerGains(mpi_rank % 2, 2, fit_results.pids_per_process);

    if (mpi_rank < 2)
    {
        save_gain_table(ensure_path_exists(gain_table_file_path(program_conf, real_radiometer)),
                        real_radiometer,
                        gTable);
        save_gain_table(ensure_path_exists(slow_gain_table_file_path(program_conf, real_radiometer)),
                        real_radiometer,
                        slowTable);
    }

    log->decrease_indent();
    log->info("Quitting module smoothGains");
}
