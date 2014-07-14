#include "mpi_processes.hpp"
#include <cstddef>

std::vector<int>
splitOdsIntoMpiProcesses(int numOfMpiProcesses,
			 const std::vector<Od_t> & list_of_ods)
{
    const int detectorsPerRadiometer = 2;
    std::vector<int> result;

    // Compute ranges
    int numOfOds = list_of_ods.size();

    for (size_t prox = 0;
	 prox < numOfMpiProcesses / detectorsPerRadiometer;
	 ++prox)
    {
	int step = detectorsPerRadiometer * numOfOds / numOfMpiProcesses + 1;
	int start = step * prox;
	int stop = start + step;
	if (stop > numOfOds)
	    stop = numOfOds;

	int nIds = 0;
	for (int intOd = start; intOd < stop; ++intOd) {
	    nIds += list_of_ods.at(intOd).num_of_pointings;
	}

	result.push_back(nIds);
    }

    return result;
}
