#include "smooth_interface.h"
#include <vector>

#include <iostream>

int main(int argc, const char ** argv)
{
    const int npids = 45000;
    int * pid = new int[npids];
    double * gain = new double[npids];
    double * dipole = new double[npids];

    for (int i=0; i<npids; i++) {
        pid[i] = i;
        gain[i] = 1.;
        dipole[i] = 1.;
    }

   int windowLenMinima = 1800;
   int windowLenMaxima = 400;
   int windowLenSlowSmoothing = 150;
   double percentSlowVariations = 0.995;
   double minRangeDipole = 3.0e-3;
   double maxRangeDipole = 3.5e-3;

    const int nJumps = 4;

    int jumps[nJumps] = {0,3267,5044,36671};
    double outputGain[npids];
    smoothGains(npids, pid, gain, dipole, outputGain, windowLenMinima,windowLenMaxima,windowLenSlowSmoothing,percentSlowVariations,minRangeDipole,maxRangeDipole, jumps, nJumps);

    std::cout << "pid" << "," << "gain" << std::endl;
    for (int i=0; i<npids; i++)
{
    std::cout << pid[i] << "," << outputGain[i] << std::endl;
}


    return 0;
}
