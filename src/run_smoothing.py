import ctypes
import numpy as np
import pandas as pd

import os
import os.path

FOLDER = os.path.dirname(__file__)
if not FOLDER:
    FOLDER = "./"

smoothlib = ctypes.CDLL(os.path.join(FOLDER, "libsmooth_interface.so"))

smoothGainsName = ["000000000001ad10 T _Z11smoothGainsiPiPdS0_S0_iiidddS_i"]
smoothGainsName = smoothGainsName[0].split()[2]

smoothGains = smoothlib[smoothGainsName]

c_float_p = ctypes.POINTER(ctypes.c_double)
c_int_p = ctypes.POINTER(ctypes.c_int64)

def smooth(raw_gains, chtag):

    dip = pd.read_csv(os.path.join(FOLDER, "dip_%s.csv" % chtag)).set_index("pointingID")

    # input arrays
    dipole = np.array(dip.deltaT.reindex(raw_gains.index).fillna(method="ffill").fillna(method="bfill")).astype(np.float64)
    pids = np.array(raw_gains.index).astype(np.int32)
    gains = np.array(raw_gains).astype(np.float64)

    # output array
    smoothed_gains = gains.copy()
    smoothed_gains[:] = 0.

    # channel parameters
    smoothing_parameters = pd.read_csv(os.path.join(FOLDER, "smoothing_parameters.csv")).set_index(["horn", "rad"])

    horn = int(chtag[3:5])
    rad = 0 if chtag.endswith("M") else 1

    windowLenMinima = int(smoothing_parameters.ix[horn, rad].smooth_window_near_dipole_min)
    windowLenMaxima  = int(smoothing_parameters.ix[horn, rad].smooth_window_near_dipole_max)
    windowLenSlowSmoothing = int(smoothing_parameters.ix[horn, rad].smooth_window_length)
    percentSlowVariations = np.float64(smoothing_parameters.ix[horn, rad].slow_var_percentile)
    minRangeDipole = np.float64(smoothing_parameters.ix[horn, rad].dipole_range_min_value)
    maxRangeDipole = np.float64(smoothing_parameters.ix[horn, rad].dipole_range_max_value)

    dx12_jumps = {
        "LFI18M": [0,3826,37290,37766],
        "LFI18S": [0,25079,43906],
        "LFI19M": [0,3267,5044,36671],
        "LFI19S": [0,5041,15011,15725,21859,26439],
        "LFI20M": [0,10992,17104,44023],
        "LFI20S": [0,16063,24854,25630],
        "LFI21M": [0,27010,43270,43928],
        "LFI21S": [0,5032,8614],
        "LFI22M": [0,3328,5040,10897,19415],
        "LFI22S": [0,3018,5025,42339,43997],
        "LFI23M": [0,2291],
        "LFI23S": [0,5055,10897,16801,43273,43813],
        "LFI24M": [0,5035,15971],
        "LFI24S": [0,5037],
        "LFI25M": [0,5037],
        "LFI25S": [0,25206],
        "LFI26M": [0,5039],
        "LFI26S": [0,5034],
        "LFI27M": [0,5044,10897],
        "LFI27S": [0,5048,10897,24764,35713],
        "LFI28M": [0,5038,10992],
        "LFI28S": [0,5039,10897],
    }
    jumps = np.array(dx12_jumps[chtag]).astype(np.int32)
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
            ctypes.c_double(maxRangeDipole),
            jumps.ctypes.data_as(c_int_p),
            len(jumps)
            )

    return pd.Series(smoothed_gains, index=raw_gains.index)
