#include <vector>
#include <iostream>
#include <numeric>

// Each array contains 9 x 2 = 18 elements
std::vector<int> pointingIds { 
    1, 2, 3, 4, 10, 11, 12, 20, 30, 31, 32, 33, 34, 40, 41, 42, 50, 51 };
std::vector<double> gain { 
    0, 0.1, 0.2, 0.3, 10, 10.1, 10.2, 20, 30, 30.1, 30.2, 30.3, 30.4, 40, 40.1, 40.2, 50, 50.1 };
std::vector<double> offset { 
    1, 1.1, 1.2, 1.3, 11, 11.1, 11.2, 21, 31, 31.1, 31.2, 31.3, 31.4, 41, 41.1, 41.2, 51, 51.1 };

////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////

void
selectRadiometerGains(int detectorIdIdx,
		      size_t detectorIdsSize,
		      const std::vector<int> & nIdsRange)
{
    // Select values
    std::vector<int> tmpPointingIds;
    std::vector<double> tmpGain;
    std::vector<double> tmpOffset;
    for (size_t idx = 0; idx < nIdsRange.size(); ++idx)
    {
	int offsetIdx = detectorIdIdx*nIdsRange.at(idx);
	int startPoint = 0;

	for (size_t intIdx = 0; intIdx < idx; ++intIdx)
	    startPoint += nIdsRange.at(intIdx)*static_cast<int>(detectorIdsSize);

	for (int intIdx = 0; intIdx < nIdsRange[idx]; ++intIdx)
	{
	    tmpPointingIds.push_back(pointingIds.at(startPoint + intIdx + offsetIdx));
	    tmpGain.push_back(gain.at(startPoint + intIdx + offsetIdx));
	    tmpOffset.push_back(offset.at(startPoint + intIdx + offsetIdx));
	}
    }

    std::cout << "detectorIdIdx = " << detectorIdIdx << std::endl;
    std::cout << "detectorIdsSize = " << detectorIdsSize << std::endl;
    print_vector(std::cout, "tmpPointingIds", tmpPointingIds);
    print_vector(std::cout, "tmpGain", tmpGain);
    print_vector(std::cout, "tmpOffset", tmpOffset);
}

////////////////////////////////////////////////////////////////////////////////

void
new_selectRadiometerGains(int detectorIdIdx,
			  size_t detectorIdsSize,
			  const std::vector<int> & nIdsRange)
{
    // Select values
    std::vector<int> tmpPointingIds;
    std::vector<double> tmpGain;
    std::vector<double> tmpOffset;

    size_t startPoint = 0;
    for (auto chunkSize : nIdsRange)
    {
	size_t chunkOffset = detectorIdIdx * chunkSize;
	size_t chunkStart = startPoint + chunkOffset;

	tmpPointingIds.insert(tmpPointingIds.end(),
			      pointingIds.begin() + chunkStart,
			      pointingIds.begin() + chunkStart + chunkSize);
	tmpGain.insert(tmpGain.end(),
		       gain.begin() + chunkStart,
		       gain.begin() + chunkStart + chunkSize);
	tmpOffset.insert(tmpOffset.end(),
			 offset.begin() + chunkStart,
			 offset.begin() + chunkStart + chunkSize);

	startPoint += detectorIdsSize * chunkSize;
    }

    std::cout << "detectorIdIdx = " << detectorIdIdx << std::endl;
    std::cout << "detectorIdsSize = " << detectorIdsSize << std::endl;
    print_vector(std::cout, "tmpPointingIds", tmpPointingIds);
    print_vector(std::cout, "tmpGain", tmpGain);
    print_vector(std::cout, "tmpOffset", tmpOffset);
}

////////////////////////////////////////////////////////////////////////////////

int main()
{
    std::vector<int> nIdsRange { 2, 4, 3 }; // Sum = 9
    selectRadiometerGains(0, 2, nIdsRange);
    selectRadiometerGains(1, 2, nIdsRange);

    new_selectRadiometerGains(0, 2, nIdsRange);
    new_selectRadiometerGains(1, 2, nIdsRange);
}
