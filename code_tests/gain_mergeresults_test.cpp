#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <mpi.h>

std::vector<int> pointingIds;
std::vector<double> gain;
std::vector<double> offset;

template<typename T>
void
print_vector(std::ostream & stream,
	     const char * name,
	     const std::vector<T> & v)
{
    stream << name << " = [ ";
    for(auto value : v) {
	stream << value << ' ';
    }
    stream << "]" << std::endl;
}

void mergeResults()
{
  // Get MPI infos
  int rankMPI = MPI::COMM_WORLD.Get_rank();
  int sizeMPI = MPI::COMM_WORLD.Get_size();

  // Set gain and pid arrays
  int length = static_cast<int>(pointingIds.size());
  int * lengths = new int[sizeMPI];
  std::fill (lengths,lengths+sizeMPI,0);
  for (int jj=0; jj<sizeMPI; ++jj)
    {
      if (jj==rankMPI)
	lengths[jj]=length;
      else
	lengths[jj]=0.0;
    }

  MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, lengths, sizeMPI, MPI::INT, MPI::SUM);
  MPI::COMM_WORLD.Barrier();

  MPI::COMM_WORLD.Allreduce(&length, &length, 1, MPI::INT, MPI::SUM);
  MPI::COMM_WORLD.Barrier();

  // Retrieve data from other processes
  int startPid=0;
  for (int jj=0; jj<rankMPI; ++jj)
    startPid+=lengths[jj];

  int * pids = new int[length];
  std::fill (pids, pids+length, 0);
  double * gainArr = new double[length];
  std::fill (gainArr, gainArr+length, 0.0);
  double * offsetArr = new double[length];
  std::fill (offsetArr, offsetArr+length, 0);

  for (unsigned int jj=0; jj<pointingIds.size(); ++jj)
    {
      pids[startPid+jj]=pointingIds[jj];
      gainArr[startPid+jj]=gain[jj];
      offsetArr[startPid+jj]=offset[jj];
    }

  MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, pids, length, MPI::INT, MPI::SUM);
  MPI::COMM_WORLD.Barrier();

  MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, gainArr, length, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Barrier();

  MPI::COMM_WORLD.Allreduce(MPI::IN_PLACE, offsetArr, length, MPI::DOUBLE, MPI::SUM);
  MPI::COMM_WORLD.Barrier();

  std::vector<int>().swap(pointingIds);
  std::vector<double>().swap(gain);
  std::vector<double>().swap(offset);

  pointingIds.assign(pids, pids+length);
  gain.assign(gainArr, gainArr+length);
  offset.assign(offsetArr, offsetArr+length);

  delete [] lengths;
  delete [] gainArr;
  delete [] pids;
  delete [] offsetArr;
}

////////////////////////////////////////////////////////////////////////////////

void improved_mergeResults()
{
  int rankMPI = MPI::COMM_WORLD.Get_rank();
  int sizeMPI = MPI::COMM_WORLD.Get_size();

  // Collect the number of pointing periods processed by each MPI
  // process into the "lengths" vector, and the total number of
  // periods into "overallLength".
  int length = static_cast<int>(pointingIds.size());
  int overallLength;
  MPI::COMM_WORLD.Allreduce(&length, &overallLength, 1, MPI::INT, MPI::SUM);

  std::vector<int> lengths(sizeMPI);
  MPI::COMM_WORLD.Gather(&length, 1, MPI::INT, 
			 lengths.data(), 1, MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(lengths.data(), lengths.size(), MPI::INT, 0);

  // Now use MPI's "gatherv" function to concatenate pIDs, gains, and
  // offsets into "overallPointings, overallGains, and overallOffsets.
  std::vector<int> displacement(sizeMPI);
  for(int i = 0; i < sizeMPI; ++i) {
      if(i == 0)
	  displacement[i] = 0;
      else
	  displacement[i] = displacement.at(i - 1) + lengths.at(i - 1);
  }

  std::vector<int> overallPointings(overallLength);
  std::vector<double> overallGains(overallLength);
  std::vector<double> overallOffsets(overallLength);

  MPI::COMM_WORLD.Gatherv(pointingIds.data(), pointingIds.size(), MPI::INT,
			  overallPointings.data(), lengths.data(),
			  displacement.data(), MPI::INT, 0);

  MPI::COMM_WORLD.Gatherv(gain.data(), gain.size(), MPI::DOUBLE,
			  overallGains.data(), lengths.data(),
			  displacement.data(), MPI::DOUBLE, 0);

  MPI::COMM_WORLD.Gatherv(offset.data(), offset.size(), MPI::DOUBLE,
			  overallOffsets.data(), lengths.data(), 
			  displacement.data(), MPI::DOUBLE, 0);

  if(rankMPI == 0) {
      pointingIds = overallPointings;
      gain = overallGains;
      offset = overallOffsets;
  } else {
      pointingIds.resize(overallLength);
      gain.resize(overallLength);
      offset.resize(overallLength);
  }

  // So far overallPointings, overallGains, and overallOffsets have
  // been set up in the root process only. Broadcast them to every
  // other process.
  MPI::COMM_WORLD.Bcast(pointingIds.data(), pointingIds.size(), MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(gain.data(), gain.size(), MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(offset.data(), offset.size(), MPI::DOUBLE, 0);
}

void
print_results(int rankMPI)
{
    std::stringstream s;
    s << "Rank #" << rankMPI << ":\n";
    print_vector(s, "\tpointingIds", pointingIds);
    print_vector(s, "\tgain", gain);
    print_vector(s, "\toffset", offset);
    std::cout << s.str() << std::endl;
}

typedef void callback_t();

int
do_test(callback_t *fn)
{
    int sizeMPI = MPI::COMM_WORLD.Get_size();
    int rankMPI = MPI::COMM_WORLD.Get_rank();

    if(sizeMPI != 6) {
	if(rankMPI == 0) {
	    std::cerr << "Error: you must use *exactly* 6 MPI processes\n";
	}

	return 1;
    }

    switch(rankMPI) {
    case 0: 
	pointingIds = { 1, 2, 3, 4 };
	gain = { 0.0, 0.1, 0.2, 0.3 };
	offset = { 1.0, 1.1, 1.2, 1.3 };
	break;
    case 1: 
	pointingIds = { 5, 6, 7 };
	gain = { 10.0, 10.1, 10.2 };
	offset = { 11.0, 11.1, 11.2 };
	break;
    case 2:
	pointingIds = { 8 };
	gain = { 20.0  };
	offset = { 21.0 };
	break;
    case 3:
	pointingIds = { 9, 10, 11, 12, 13 };
	gain = { 30.0, 30.1, 30.2, 30.3, 30.4 };
	offset = { 31.0, 31.1, 31.2, 31.3, 31.4 };
	break;
    case 4:
	pointingIds = { 14, 15, 16 };
	gain = { 40.0, 40.1, 40.2 };
	offset = { 41.0, 41.1, 41.2 };
	break;
    case 5:
	pointingIds = { 17, 18, 19 };
	gain = { 50.0, 50.1, 50.2 };
	offset = { 51.0, 51.1, 51.2 };
	break;
    }

    fn();

    // Implement a round-robin strategy to make all the MPI processes
    // write their results on the screen.
    if(rankMPI == 0) {
	std::cout << "Results:\n\n";
	print_results(rankMPI);
	MPI::COMM_WORLD.Send(&rankMPI, 1, MPI::INT, rankMPI + 1, 0);
    } else {
	int dummy;
	MPI::COMM_WORLD.Recv(&dummy, 1, MPI::INT, rankMPI - 1, 0);
	print_results(rankMPI);
	if(rankMPI < sizeMPI - 1)
	    MPI::COMM_WORLD.Send(&rankMPI, 1, MPI::INT, rankMPI + 1, 0);
    }

    return 0;
}

int
main()
{
    MPI::Init(); 

    int result = do_test(mergeResults);

    MPI::COMM_WORLD.Barrier();

    if(result == 0) {
	result = do_test(improved_mergeResults);
    }

    MPI::Finalize();
    return result;
}
