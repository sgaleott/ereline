#include <iostream>
#include "gain_table.hpp"
#include "smooth_interface.h"

int smoothGains(int npids,
                int * pid,
                double * gain,
                double * dipole,
                double * outputGain,
                int windowLenMinima,
                int windowLenMaxima,
                int windowLenSlowSmoothing,
                double percentSlowVariations,
                double minRangeDipole,
                double maxRangeDipole,
                int * jumps,
                int nJumps)
{
// horn|rad|smooth_window_near_dipole_min|smooth_window_near_dipole_max|smooth_window_length|fast_variations_window_length|slow_var_percentile|fast_variations_percentile|dipole_range_min_value|dipole_range_max_value
// 19|1|1800|400|150|300|0.995|0.95|0.003|0.0035

//    int windowLenMinima = 1800;
//    int windowLenMaxima = 400;
//    int windowLenSlowSmoothing = 150;
//    double percentSlowVariations = 0.995;
//    double minRangeDipole = 3.0e-3;
//    double maxRangeDipole = 3.5e-3;
    std::vector<double> dipoleVec;
    std::vector<int> jumpsVec(jumps, jumps + nJumps);

    std::cout << pid[0] << std::endl;
    std::cout << pid[1] << std::endl;
    std::cout << gain[0] << std::endl;
    std::cout << gain[1] << std::endl;
    std::cout << dipole[0] << std::endl;
    std::cout << dipole[1] << std::endl;

    Gain_table_t outRawTable;
    for (int i = 0; i < npids; i++)
    {
        outRawTable.append({pid[i], gain[i], 0. });
        dipoleVec.push_back(dipole[i]);
    }
    std::vector<double> outputGainVec = outRawTable.gainSmoothing(windowLenMinima, windowLenMaxima,
            minRangeDipole, maxRangeDipole, jumpsVec, dipoleVec);
    for (int i = 0; i < npids; i++)
    {
        outputGain[i] = outputGainVec[i];
    }
    std::cout << outputGain[0] << std::endl;
    std::cout << outputGain[1] << std::endl;
    return 0;
}
