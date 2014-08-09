#include <mpi.h>

#include "data_binning_results.hpp"

#include "configuration.hpp"
#include "fits_object.hpp"
#include "io.hpp"
#include "logging.hpp"

#include <algorithm>

extern "C" {
#include "Healpix_3.11/chealpix.h"
};

const size_t MAX_PATH_LENGTH = 1024;

////////////////////////////////////////////////////////////////////////////////

static std::string
get_base_path(const Configuration & program_conf)
{
    return program_conf.getWithSubst("common.base_output_dir") 
        + "/binned_data/";
}

////////////////////////////////////////////////////////////////////////////////

static std::string
get_index_file_name(const Configuration & program_conf,
                    const Lfi_radiometer_t & radiometer)
{
    return ((boost::format("%s/%s_binned_data_index.fits")
             % get_base_path(program_conf)
             % radiometer.shortName())
            .str());
}

////////////////////////////////////////////////////////////////////////////////

static std::string
get_binned_data_file_name(const Configuration & program_conf,
                          const Lfi_radiometer_t & radiometer,
                          int pointing_id,
                          int od)
{
    // Example: "basedir/LFI27M/OD0091/binned_pid_od0091_000042.fits"
    return ((boost::format("%1%/%2%/OD%3$04d/binned_pid_od%3$04d_%4$06d.fits")
             % get_base_path(program_conf)
             % radiometer.shortName()
             % od
             % pointing_id)
            .str());
}

////////////////////////////////////////////////////////////////////////////////

static std::string
get_mask_file_name(const Configuration & program_conf,
                   const Lfi_radiometer_t & radiometer)
{
    return ((boost::format("%s/%s_binning_mask.fits.gz")
             % get_base_path(program_conf)
             % radiometer.shortName())
            .str());
}

////////////////////////////////////////////////////////////////////////////////

static void
save_binned_pid(const Binned_data_t & cur_binned_pid,
                const std::string & file_name)
{
    FitsObject pid_file;
    pid_file.create(ensure_path_exists(file_name));

    const std::vector<fitscolumn> pid_columns {
        { "PIXIDX", "", fitsTypeC<long>(), 1 },
        { "HITS", "", fitsTypeC<long>(), 1 },
        { "DATASUM", "", fitsTypeC<double>(), 1 },
        { "MODELAVG", "", fitsTypeC<float>(), 1 }
    };
    pid_file.insertBINtable(pid_columns, "BINDATA");
    pid_file.writeColumn(1, cur_binned_pid.pix_index);
    pid_file.writeColumn(2, cur_binned_pid.pix_num_of_hits);
    pid_file.writeColumn(3, cur_binned_pid.pix_data_sum);
    pid_file.writeColumn(4, cur_binned_pid.pix_model_mean);

    pid_file.setKey("QUALFLAG", cur_binned_pid.quality_flag,
                    "Quality flag used for masking bad data");
    pid_file.setKey("NSIDE", cur_binned_pid.nside,
                    "NSIDE used in the Healpix pixelisation");
    pid_file.setKey("PID", cur_binned_pid.pointing_id);
    pid_file.setKey("OD", cur_binned_pid.od);
    pid_file.setKey("MAXDIPOL", cur_binned_pid.max_dipole);
    pid_file.setKey("MINDIPOL", cur_binned_pid.min_dipole);
}

////////////////////////////////////////////////////////////////////////////////

static void
load_binned_pid(const std::string & file_name,
                Binned_data_t & binned_pid)
{
    FitsObject pid_file;
    pid_file.openTable(file_name);

    pid_file.getColumn(1, binned_pid.pix_index);
    pid_file.getColumn(2, binned_pid.pix_num_of_hits);
    pid_file.getColumn(3, binned_pid.pix_data_sum);
    pid_file.getColumn(4, binned_pid.pix_model_mean);

    pid_file.getKey("QUALFLAG", binned_pid.quality_flag);
    pid_file.getKey("NSIDE", binned_pid.nside);
    pid_file.getKey("PID", binned_pid.pointing_id);
    pid_file.getKey("OD", binned_pid.od);
    pid_file.getKey("MAXDIPOL", binned_pid.max_dipole);
    pid_file.getKey("MINDIPOL", binned_pid.min_dipole);
}

////////////////////////////////////////////////////////////////////////////////

static void
save_index_file(const Lfi_radiometer_t & radiometer,
                const std::vector<int> & pids,
                const std::vector<int> & ods,
                const std::vector<int> & pids_per_process,
                const std::string & index_file_path)
{
    FitsObject index_file;
    index_file.create(ensure_path_exists(index_file_path));

    const std::vector<fitscolumn> index_columns {
        { "PID", "", fitsTypeC<int>(), 1 },
        { "OD", "", fitsTypeC<int>(), 1 }
    };
    index_file.insertBINtable(index_columns, "PIDIDX");
    index_file.writeColumn(1, pids);
    index_file.writeColumn(2, ods);

    const std::vector<fitscolumn> mpi_columns {
        { "PIDSPROC", "", fitsTypeC<int>(), 1 }
    };
    index_file.insertBINtable(mpi_columns, "MPIINFO");
    index_file.writeColumn(1, pids_per_process);

    const int mpi_size = MPI::COMM_WORLD.Get_size();
    index_file.setKey("MPISIZE", mpi_size, "Number of MPI processes used");
}

////////////////////////////////////////////////////////////////////////////////

static void
load_index_file(const std::string & index_file_path,
                Lfi_radiometer_t & radiometer,
                std::vector<int> & pids,
                std::vector<int> & ods,
                std::vector<int> & pids_per_process)
{
    FitsObject index_file;
    index_file.openTable(index_file_path);

    index_file.getColumn(1, pids);
    index_file.getColumn(2, ods);

    index_file.gotoHDU(3);

    int mpi_size_from_file;
    Logger * log = Logger::get_instance();
    index_file.getKey("MPISIZE", mpi_size_from_file);

    const int mpi_size = MPI::COMM_WORLD.Get_size();
    if(mpi_size != mpi_size_from_file) {
        throw std::runtime_error(
            (boost::format("the binned data have been created using %1% "
                           "MPI processes, but now %2% processes are running")
             % mpi_size_from_file
             % mpi_size).str());
    }

    index_file.getColumn(1, pids_per_process);
}

////////////////////////////////////////////////////////////////////////////////

void Data_binning_results_t::save_to_disk(Configuration & program_conf) const
{
    Logger * log = Logger::get_instance();
    log->debug(boost::format("Entering function %1%")
               % __PRETTY_FUNCTION__);
    log->increase_indent();

    // Step 1: save all the binned data in separate files
    const std::string base_path(get_base_path(program_conf));
    
    log->debug("Saving all binned data in FITS files");
    log->increase_indent();

    std::vector<int> pids;
    std::vector<int> ods;

    for(auto const & cur_binned_pid : binned_pids) {
        const std::string file_name(
            get_binned_data_file_name(program_conf, radiometer,
                                      cur_binned_pid.pointing_id,
                                      cur_binned_pid.od));

        try {
            log->debug(boost::format("Saving pID %1% (OD %2%) in %3%")
                       % cur_binned_pid.pointing_id
                       % cur_binned_pid.od
                       % file_name);
            save_binned_pid(cur_binned_pid, file_name);
        }
        catch(std::exception & exc) {
            log->warning(boost::format("Unable to save pID %1% (%2%), "
                                       "skipping...")
                         % cur_binned_pid.pointing_id
                         % exc.what());
            continue;
        }

        pids.push_back(cur_binned_pid.pointing_id);
        ods.push_back(cur_binned_pid.od);
    }
    log->decrease_indent();
    log->debug("Binned data saved");
    
    // Step 2: save the "index file" (which contains a list of all the
    // PIDs included in this result, as well as the information about
    // the way MPI processes have split the data among them)
    const std::string index_file_name(get_index_file_name(program_conf,
                                                          radiometer));
    log->debug(boost::format("Creating index file %1% of binned data files")
               % index_file_name);
    save_index_file(radiometer, pids, ods, pids_per_process, index_file_name);

    // Step 3: save the mask
    const std::string mask_file_name(
        get_mask_file_name(program_conf, radiometer));
    log->debug(boost::format("Creating mask file %1%")
               % mask_file_name);
    write_healpix_map(mask.pixels.data(), 
                      mask.nside,
                      ("!" + mask_file_name).c_str(),
                      mask.ordering == Healpix::Ordering_t::NEST,
                      "E");

    log->decrease_indent();
    log->debug(boost::format("Quitting function %1%")
               % __PRETTY_FUNCTION__);
}

////////////////////////////////////////////////////////////////////////////////

void Data_binning_results_t::load_from_disk(const Lfi_radiometer_t &a_radiometer,
                                            Configuration & program_conf)
{
    Logger * log = Logger::get_instance();
    log->debug(boost::format("Entering function %1%")
               % __PRETTY_FUNCTION__);
    log->increase_indent();

    radiometer = a_radiometer;

    // Step 1: load the "index file"
    const std::string index_file_name(get_index_file_name(program_conf,
                                                          radiometer));
    log->debug(boost::format("Loading index file %1%")
               % index_file_name);

    std::vector<int> pids;
    std::vector<int> ods;
    load_index_file(index_file_name, radiometer,
                    pids, ods, pids_per_process);

    // Step 2: load all the FITS files containing binned data
    const std::string base_path(get_base_path(program_conf));
    
    log->debug("Loading all binned data from FITS files");
    log->increase_indent();
    for(size_t idx = 0; idx < pids.size(); ++idx) {
        const std::string file_name(
            get_binned_data_file_name(program_conf, radiometer,
                                      pids[idx], ods[idx]));

        Binned_data_t cur_binned_pid;
        try {
            log->debug(boost::format("Loading pID %1% (OD %2%) from %3%")
                       % cur_binned_pid.pointing_id
                       % cur_binned_pid.od
                       % file_name);
            load_binned_pid(file_name, cur_binned_pid);
        }
        catch(std::exception & exc) {
            log->warning(boost::format("Unable to load pID %1% (%2%), "
                                       "skipping...")
                         % cur_binned_pid.pointing_id
                         % exc.what());
            continue;
        }

        binned_pids.push_back(cur_binned_pid);
    }
    log->decrease_indent();
    log->debug("Binned data loaded");
    
    // Step 3: load the mask
    const std::string mask_file_name(
        get_mask_file_name(program_conf, radiometer));
    long nside;
    char coordsys[9];
    char ordering[9];
    const float * pixels = 
        read_healpix_map(mask_file_name.c_str(),
                         &nside, coordsys, ordering);
    const size_t num_of_pixels = 12 * nside * nside;

    mask.nside = nside;
    mask.ordering = 
        (ordering[0] == 'N') ? 
        Healpix::Ordering_t::NEST : 
        Healpix::Ordering_t::RING;
    mask.pixels.insert(mask.pixels.begin(),
                       pixels, pixels + num_of_pixels);

    log->info(boost::format("Map %1% loaded, coordsys=\"%2%\", "
                            "ordering=\"%3%\", nside=%4%")
              % mask_file_name
              % coordsys
              % ordering
              % nside);

    log->decrease_indent();
    log->debug(boost::format("Quitting function %1%")
               % __PRETTY_FUNCTION__);
}
