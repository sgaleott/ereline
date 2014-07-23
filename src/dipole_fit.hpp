#ifndef _DIPOLE_FIT_
#define _DIPOLE_FIT_

#include <vector>

#include "ahf_info.hpp"
#include "gain_table.hpp"
#include "healpix_map.hpp"
#include "planck_velocity.hpp"

struct Lfi_radiometer_t;

/* class calib iter*/
struct Dipole_fit_t
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

    Dipole_fit_t(uint32_t a_qualityFlag, int a_nSide, int a_pointingID);

    bool binData(const std::vector<double> & data,
                 const std::vector<uint32_t>& flag,
                 const std::vector<double> & theta,
                 const std::vector<double> & phi,
                 const std::vector<double> & dipole,
                 const Range_t<size_t> & index_range,
                 const std::vector<double> & sidelobes);
    bool fitData(const std::vector<float> & maskMap);
    bool fit(const std::vector<double> & data,
             const std::vector<uint32_t> & flag,
             const std::vector<double> & theta,
             const std::vector<double> & phi,
             const std::vector<double> & dipole,
             const Range_t<size_t> & index_range,
             const std::vector<float> & maskMap,
             const std::vector<double> & sidelobes);

    void setPixSumDipole(const std::vector<float> & inpArr);

    double getDipoleVariance() const;
    double getMaxDipole() const;
    double getMinDipole() const;

    void unload();
};

class Configuration;
struct Sqlite_connection_t;
struct Dipole_fit_results_t;

void run_dipole_fit(Sqlite_connection_t & ucds,
                    const Lfi_radiometer_t & rad,
                    Configuration & program_conf,
                    Configuration & storage_conf,
                    const std::vector<Pointing_t> & list_of_pointings,
                    Dipole_fit_results_t & result);

#endif
