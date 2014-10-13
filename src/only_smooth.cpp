#include "smooth_interface.h"
#include <vector>

#include <iostream>
#include "LFI19M_pids_gains_dipole.cpp"

int main(int argc, const char ** argv)
{
    
    double outputGain[npids];
    smoothGains(npids, pid, gain, dipole, outputGain);

    std::cout << "pid" << "," << "gain" << std::endl;
    for (int i=0; i<npids; i++)
{
    std::cout << pid[i] << "," << outputGain[i] << std::endl;
}


    return 0;
}
