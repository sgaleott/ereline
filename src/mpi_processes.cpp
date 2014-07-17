#include "mpi_processes.hpp"
#include <mpi.h>
#include <cstddef>

void
splitOdsIntoMpiProcesses(int numOfMpiProcesses,
			 const std::vector<Od_t> & list_of_ods,
			 std::vector<Data_range_t> & data_range)
{
    const int detectorsPerRadiometer = 2;

    // Compute ranges
    const int numOfOds = list_of_ods.size();

    data_range.clear();
    for (size_t prox = 0;
	 prox < numOfMpiProcesses / detectorsPerRadiometer;
	 ++prox)
    {
	int step = detectorsPerRadiometer * numOfOds / numOfMpiProcesses + 1;
	int start = step * prox;
	int stop = start + step;
	if (stop > numOfOds)
	    stop = numOfOds;

	Data_range_t cur_range;

	if(start < numOfOds && stop > start) {
	    cur_range = { 
		{ list_of_ods.at(start).od, list_of_ods.at(stop - 1).od },
		{ list_of_ods.at(start).first_pointing_id, list_of_ods.at(stop - 1).last_pointing_id },
		stop - start,
		0
	    };
	} else {
	    cur_range = { { -1, -1 }, { -1, -1 }, 0, 0 };
	}

	for (int od_idx = start; od_idx < stop; ++od_idx) {
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
