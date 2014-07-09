#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <cassert>

const std::vector<int> pids_per_od { 10, 14, 9, 11, 11, 8, 13, 17, 12 };

void
ahfInfos_selectProcessRange(int sizeMPI, size_t detectorIdsSize,
			    std::vector<int> & nIdsRange)
{
    assert(sizeMPI % detectorIdsSize == 0);

    // Compute ranges
    int odNumber = pids_per_od.size();
    const int step = detectorIdsSize*odNumber/sizeMPI + 1;
    for (size_t prox=0; prox<sizeMPI / detectorIdsSize; ++prox)
    {
	int start = step * prox;
	int stop = step * (prox + 1);
	if (stop > odNumber)
	    stop = odNumber;

	int nIds = 0;
	for (int intOd = start; intOd < stop; ++intOd)
	    nIds += pids_per_od.at(intOd);

	nIdsRange.push_back(nIds);
    }
}

template<typename T>
void
print_vector(const char * name, const std::vector<T> & v)
{
    std::cout << name << " = [ ";
    for(auto value : v) {
	std::cout << value << ' ';
    }
    std::cout << "]" << std::endl;
}

int main(int argc, const char ** argv)
{
    std::vector<int> nIdsRange;
    int sizeMPI;
    const int detectorIdsSize = 2;

    if(argc != 2) {
	std::cerr << "Usage: select_proc_range NUM_OF_MPI_PROCESSES\n";
	return 1;
    }
    sizeMPI = std::atoi(argv[1]);

    std::cout << "sizeMPI = " << sizeMPI << std::endl;
    std::cout << "detectorIdsSize = " << detectorIdsSize << std::endl;
    print_vector("pids_per_od", pids_per_od);

    ahfInfos_selectProcessRange(sizeMPI,
				detectorIdsSize,
				nIdsRange);
    print_vector("nIdsRange (result of ahfInfos::selectProcessRange)",
		 nIdsRange);

    return 0;
}
