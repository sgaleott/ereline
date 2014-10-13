#ifndef SMOOTH_INTERFACE_HPP
#define SMOOTH_INTERFACE_HPP

int smoothGains(int npids, int * pid, double * gain, double * dipole, double * outputGain,
                int windowLenMinima,
                int windowLenMaxima,
                int windowLenSlowSmoothing,
                double percentSlowVariations,
                double minRangeDipole,
                double maxRangeDipole);

#endif
