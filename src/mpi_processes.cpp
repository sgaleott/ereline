#include <mpi.h>

#include "mpi_processes.hpp"
#include "configuration.hpp"
#include "logging.hpp"
#include <cstddef>

void
splitOdsIntoMpiProcesses(int numOfMpiProcesses,
                         const std::vector<Od_t> & list_of_ods,
                         std::vector<Data_range_t> & data_range)
{
    const int detectorsPerRadiometer = 2;
    const int num_of_chunks = numOfMpiProcesses / detectorsPerRadiometer;
    // Compute ranges
    const int numOfOds = list_of_ods.size();

    data_range.clear();
    for (size_t prox = 0; prox < num_of_chunks; ++prox)
    {
        const int start = int(std::round(float(prox * numOfOds) / num_of_chunks));
        const int stop = int(std::round(float((prox + 1) * numOfOds) / num_of_chunks)) - 1;

        Data_range_t cur_range;

        if(start < numOfOds && stop >= start) {
            cur_range = {
                { list_of_ods.at(start).od, list_of_ods.at(stop).od },
                { list_of_ods.at(start).first_pointing_id, list_of_ods.at(stop).last_pointing_id },
                stop - start + 1,
                0
            };
        } else {
            cur_range = { { -1, -1 }, { -1, -1 }, 0, 0 };
        }

        for (int od_idx = start; od_idx <= stop; ++od_idx) {
            cur_range.num_of_pids += list_of_ods.at(od_idx).num_of_pointings;
        }

        data_range.push_back(cur_range);
    }
}

////////////////////////////////////////////////////////////////////////////////

Range_t<int>
range_of_ods_for_this_process(const std::vector<Od_t> & list_of_ods)
{
    const int rank = MPI::COMM_WORLD.Get_rank();
    const int size = MPI::COMM_WORLD.Get_size();

    std::vector<Data_range_t> pids_per_process;
    splitOdsIntoMpiProcesses(size, list_of_ods, pids_per_process);

    return Range_t<int> { pids_per_process.at(rank).od_range.start,
            pids_per_process.at(rank).od_range.end };
}

////////////////////////////////////////////////////////////////////////////////

/* Assuming that each MPI Process has its own local chunk of the three vectors
 * "pointingIds", "vec1", and "vec2", this function gathers them into three
 * large vectors. The vector "pointingIds" contains the number of pIDs
 * that have been handled by each MPI process, and it is used to determine
 * how many elements must be taken from "vec1" and "vec2". */
void
merge_tables(std::vector<int> & pointingIds,
             std::vector<double> & vec1,
             std::vector<double> & vec2)
{
    int rankMPI = MPI::COMM_WORLD.Get_rank();
    int sizeMPI = MPI::COMM_WORLD.Get_size();

    // Collect the number of pointing periods processed by each MPI
    // process into the "lengths" vector, and the total number of
    // periods into "overallLength".
    const int length = pointingIds.size();
    int overallLength;
    MPI::COMM_WORLD.Allreduce(&length, &overallLength, 1, MPI::INT, MPI::SUM);

    std::vector<int> lengths(sizeMPI);
    MPI::COMM_WORLD.Gather(&length, 1, MPI::INT, lengths.data(), 1, MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(lengths.data(), lengths.size(), MPI::INT, 0);

    // Now use MPI's "gatherv" function to concatenate pIDs, gains, and
    // offsets into "overallPointings, overallGains, and overallOffsets.
    std::vector<int> displacement(sizeMPI);
    displacement[0] = 0;
    for(int i = 1; i < sizeMPI; ++i) {
        displacement[i] = displacement.at(i - 1) + lengths.at(i - 1);
    }

    std::vector<int> overallPointings(overallLength);
    std::vector<double> overall1(overallLength);
    std::vector<double> overall2(overallLength);

    MPI::COMM_WORLD.Gatherv(pointingIds.data(), pointingIds.size(), MPI::INT,
                            overallPointings.data(), lengths.data(),
                            displacement.data(), MPI::INT, 0);

    MPI::COMM_WORLD.Gatherv(vec1.data(), vec1.size(), MPI::DOUBLE,
                            overall1.data(), lengths.data(),
                            displacement.data(), MPI::DOUBLE, 0);

    MPI::COMM_WORLD.Gatherv(vec2.data(), vec2.size(), MPI::DOUBLE,
                            overall2.data(), lengths.data(),
                            displacement.data(), MPI::DOUBLE, 0);

    if(rankMPI == 0) {
        pointingIds = overallPointings;
        vec1 = overall1;
        vec2 = overall2;
    } else {
        pointingIds.resize(overallLength);
        vec1.resize(overallLength);
        vec2.resize(overallLength);
    }

    // So far overallPointings, overallGains, and overallOffsets have
    // been set up in the root process only. Broadcast them to every
    // other process.
    MPI::COMM_WORLD.Bcast(pointingIds.data(), pointingIds.size(), MPI::INT, 0);
    MPI::COMM_WORLD.Bcast(vec1.data(), vec1.size(), MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(vec2.data(), vec2.size(), MPI::DOUBLE, 0);
}

////////////////////////////////////////////////////////////////////////////////

void
get_local_data_range(int mpi_rank,
                     int mpi_size,
                     const std::vector<Pointing_t> & list_of_pointings,
                     std::vector<int> & num_of_pids,
                     Data_range_t & data_range)
{
    Logger * log = Logger::get_instance();

    auto list_of_ods = build_od_list(list_of_pointings);

    std::vector<Data_range_t> list_of_data_ranges;
    splitOdsIntoMpiProcesses(mpi_size, list_of_ods, list_of_data_ranges);
    log->info(boost::format("The data to analyze will be split in %1% chunks "
                            "(there are %2% MPI jobs running)")
              % list_of_data_ranges.size()
              % mpi_size);

    num_of_pids.resize(list_of_data_ranges.size());
    for(size_t idx = 0; idx < list_of_data_ranges.size(); ++idx) {
        num_of_pids[idx] = list_of_data_ranges[idx].num_of_pids;
    }

    data_range = list_of_data_ranges.at(mpi_rank / 2);
    log->info(boost::format("Data range for this MPI process: ODs [%1%, %2%] "
                            "(pointings [%3%, %4%]), number of pointings: %5%")
              % data_range.od_range.start
              % data_range.od_range.end
              % data_range.pid_range.start
              % data_range.pid_range.end
              % data_range.num_of_pids);
}

////////////////////////////////////////////////////////////////////////////////

Lfi_radiometer_t
radiometer_to_use(int mpi_rank,
                  const Lfi_radiometer_t & user_rad,
                  Configuration & program_conf,
                  Configuration & storage_conf)
{
    Logger * log = Logger::get_instance();

    /* The way MPI processes are split during binning is the
     * following:
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
