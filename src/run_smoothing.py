import ctypes
import numpy as np
import pandas as pd

import os
FOLDER = os.path.dirname(__file__)

smoothlib = ctypes.CDLL("./libsmooth_interface.so")

smoothGainsName = get_ipython().getoutput(u'nm -D libsmooth_interface.so | grep smoothGains')
smoothGainsName = smoothGainsName[0].split()[2]

smoothGains = smoothlib[smoothGainsName]

c_float_p = ctypes.POINTER(ctypes.c_double)
c_int_p = ctypes.POINTER(ctypes.c_int64)

import pyfits

def smooth(raw_gains, chtag):

    dip = pd.read_csv(FOLDER + "dip_%s.csv" % chtag).set_index("pointingID")

    # input arrays
    dipole = np.array(dip.deltaT.reindex(raw_gains.index).fillna(method="ffill").fillna(method="bfill")).astype(np.float64)
    pids = np.array(raw_gains.index).astype(np.int32)
    gains = np.array(raw_gains).astype(np.float64)

    # output array
    smoothed_gains = gains.copy()
    smoothed_gains[:] = 0.

    # channel parameters
    smoothing_parameters = pd.read_csv(FOLDER + "smoothing_parameters.csv").set_index(["horn", "rad"])

    horn = int(chtag[3:5])
    rad = 0 if chtag.endswith("M") else 1

    windowLenMinima = int(smoothing_parameters.ix[horn, rad].smooth_window_near_dipole_min)
    windowLenMaxima  = int(smoothing_parameters.ix[horn, rad].smooth_window_near_dipole_max)
    windowLenSlowSmoothing = int(smoothing_parameters.ix[horn, rad].smooth_window_length)
    percentSlowVariations = np.float64(smoothing_parameters.ix[horn, rad].slow_var_percentile)
    minRangeDipole = np.float64(smoothing_parameters.ix[horn, rad].dipole_range_min_value)
    maxRangeDipole = np.float64(smoothing_parameters.ix[horn, rad].dipole_range_max_value)

    print maxRangeDipole

    smoothGains(len(raw_gains),
            pids.ctypes.data_as(c_int_p),
            gains.ctypes.data_as(c_float_p),
            dipole.ctypes.data_as(c_float_p),
            smoothed_gains.ctypes.data_as(c_float_p),
            windowLenMinima,
            windowLenMaxima,
            windowLenSlowSmoothing,
            ctypes.c_double(percentSlowVariations),
            ctypes.c_double(minRangeDipole),
            ctypes.c_double(maxRangeDipole))

    return pd.Series(smoothed_gains, index=raw_gains.index)
