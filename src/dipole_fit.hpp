#ifndef _DIPOLE_FIT_
#define _DIPOLE_FIT_

#include <vector>

#include "ahf_info.hpp"
#include "planck_velocity.hpp"

struct LfiRadiometer;

/* class calib iter*/
struct dipoleFit
{
    int qualityFlag;
    int nSide;
    int pointingID;

    double gainv;
    double offset;

    double maxDipole;
    double minDipole;

    std::vector<int> pixIndex;
    std::vector<int> pixSumHits;
    std::vector<double> pixSumData;
    std::vector<float> pixSumDipole;
    std::vector<float> inputMap;

    dipoleFit(uint32_t a_qualityFlag, int a_nSide, int a_pointingID);

    ~dipoleFit(){};

    bool binData(const std::vector<double> & data, 
		 const std::vector<uint32_t>& flag,
		 const std::vector<double> & theta, 
		 const std::vector<double> & phi, 
		 const std::vector<double> & dipole,
		 const std::vector<size_t> & pidRange, 
		 const std::vector<double> & sidelobes);
    bool fitData(const std::vector<float> & maskMap);
    bool fit(const std::vector<double> & data, 
	     const std::vector<uint32_t> & flag,
	     const std::vector<double> & theta, 
	     const std::vector<double> & phi, 
	     const std::vector<double> & dipole,
	     const std::vector<size_t> & pidRange, 
	     const std::vector<float> & maskMap, 
	     const std::vector<double> & sidelobes);

    void setGainV(double a_gainv);
    void setOffset(double a_offset);
    void setPixSumDipole(const std::vector<float> & inpArr);

    double getGainV() const;
    double getGain() const;
    double getOffset() const;

    int getPointingID() const;
    int getNSide() const;

    const std::vector<int> & getPixIndex() const;
    const std::vector<double> & getPixSumData() const;
    const std::vector<int> & getPixSumHits() const;
    const std::vector<float> & getPixSumDipole() const;

    double getDipoleVariance() const;
    double getMaxDipole() const;
    double getMinDipole() const;

    void unload();
};

class Configuration;

void run_dipole_fit(const LfiRadiometer & rad,
		    Configuration & program_conf,
		    Configuration & storage_conf,
		    const std::vector<Pointing_t> & list_of_pointings);

#endif
