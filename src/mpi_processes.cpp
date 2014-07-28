#include "mpi_processes.hpp"
#include <mpi.h>
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
