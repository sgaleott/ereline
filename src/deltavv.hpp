#ifndef _DELTAVONV_
#define _DELTAVONV_

#include <vector>

#include "dipole_fit.hpp"

struct Delta_vv_state_t {
    int pid;
    double gain;
    double dipole;
};

/* class calib iter*/
struct Delta_vv_t {
    std::vector<int> pid;
    std::vector<double> gain;
    std::vector<double> dipole;

    Delta_vv_t();

    void hybridFit(const std::vector<double> & subData,
                   const std::vector<double> & subSensor);

    std::vector<double> eval(const std::vector<int> & pointingIDs) const;

    void appendPidValue (int a_pid);
    void appendGainValue (double a_gain);
    void appendDipoleValue (double a_dipole);

    void append(const Delta_vv_state_t & state) {
        appendPidValue(state.pid);
        appendGainValue(state.gain);
        appendDipoleValue(state.dipole);
    }

    void mergeResults();
    void selectRadiometerGains(int detectorIdIdx,
                               size_t detectorIdsSize,
                               const std::vector<int> & nIdsRange);
};

#endif
