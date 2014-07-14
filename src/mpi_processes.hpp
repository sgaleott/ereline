#ifndef MPI_PROCESSES_HPP
#define MPI_PROCESSES_HPP

#include "ahf_info.hpp"

std::vector<int>
splitOdsIntoMpiProcesses(int numOfMpiProcesses,
			 const std::vector<Od_t> & list_of_ods);

#endif

