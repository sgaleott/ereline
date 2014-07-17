#ifndef MPI_PROCESSES_HPP
#define MPI_PROCESSES_HPP

#include "ahf_info.hpp"
#include "datatypes.hpp"

struct Data_range_t {
    Range_t<int> od_range;
    Range_t<int> pid_range;

    int num_of_ods; // This might be smaller than (od_range.end - od_range.start + 1) !
    int num_of_pids;

    // This is useful for unit tests
    bool operator!=(const Data_range_t & other) const {
	return (od_range != other.od_range) 
	    || (pid_range != other.pid_range) 
	    || num_of_ods != other.num_of_ods 
	    || num_of_pids != other.num_of_pids;
    }
};

inline std::ostream & 
operator<<(std::ostream & os, Data_range_t const & yt)
{
    os << "[ OD range: " << yt.od_range 
       << ", pID range: " << yt.pid_range 
       << ", #OD: " << yt.num_of_ods
       << ", #pID: " << yt.num_of_pids
       << "]";

    return os;
}

void
splitOdsIntoMpiProcesses(int numOfMpiProcesses,
			 const std::vector<Od_t> & list_of_ods,
			 std::vector<Data_range_t> & data_range);

Range_t<int>
range_of_ods_for_this_process(const std::vector<Od_t> & list_of_ods);

#endif

