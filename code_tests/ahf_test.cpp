#include <vector>
#include <sstream>
#include <iostream>
#include <mpi.h>

std::vector<int> nIdsRange;
const std::vector<int> pids_per_od { 10, 14, 9, 11, 11, 8, 13, 17, 12 };

void
selectProcessRange(int sizeMPI, size_t detectorIdsSize)
{
    // Compute ranges
    int odNumber = pids_per_od.size();
    for (size_t prox=0; prox<sizeMPI/detectorIdsSize; ++prox)
    {
	int step = detectorIdsSize * odNumber / sizeMPI + 1;
	int start = step * prox;
	int stop = start + step;
	if (stop > odNumber)
	    stop = odNumber;

	int nIds = 0;
	for (int intOd=start; intOd<stop; ++intOd)
	{
	    nIds += pids_per_od.at(intOd);
	}
	nIdsRange.push_back(nIds);
    }
}

std::vector<int> 
copyProcessRange(int detectorIdsSize)
{
    // Get MPI infos
    int rankMPI = MPI::COMM_WORLD.Get_rank();
    int sizeMPI = MPI::COMM_WORLD.Get_size();

    // Set zero array
    int * nIds = new int[sizeMPI/detectorIdsSize];
    if (rankMPI != 0)
	for (int idx=0; idx<sizeMPI/detectorIdsSize; ++idx)
	{
	    nIds[idx]=0;
	}
    else
    {
	for (size_t idx=0; idx<nIdsRange.size(); ++idx)
	    nIds[idx]=nIdsRange[idx];
    }

    // Retrieve data from other processes
    MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, nIds, sizeMPI/detectorIdsSize, MPI::INT, MPI::SUM);
    MPI::COMM_WORLD.Barrier();

    // Set final array
    if (rankMPI != 0)
	for (int idx=0; idx<sizeMPI/detectorIdsSize; ++idx)
	    nIdsRange.push_back(nIds[idx]);

    return nIdsRange;
}

template<typename T>
void
print_vector(const char * name, const std::vector<T> & v)
{
    std::stringstream s;
    int rankMPI = MPI::COMM_WORLD.Get_rank();
    s << "Rank #" << rankMPI << ": ";
    s << name << " = [ ";
    for(auto value : v) {
	s << value << ' ';
    }
    s << "]" << std::endl;

    std::cout << s.str();
}

int main()
{
    MPI::Init(); 
    int sizeMPI = MPI::COMM_WORLD.Get_size();
    int rankMPI = MPI::COMM_WORLD.Get_rank();

    if (rankMPI == 0)
	selectProcessRange(sizeMPI, 2);

    MPI::COMM_WORLD.Barrier(); 
    std::vector<int> idsRange = copyProcessRange(2);

    print_vector("idsRange", idsRange);

    MPI::Finalize();
}
