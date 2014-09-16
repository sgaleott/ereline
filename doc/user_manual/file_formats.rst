.. Description of the file formats used by Ereline

File formats used by Ereline
============================

In this section we provide a description of the format of the files
used by Ereline.

Pointing information
--------------------

Pointing files can be provided in two formats: FITS files and
*Squeezer* files. The program is able to automatically detect the kind
of file. The code assumes that each files contains data for one OD.

FITS files need to have the following columns:

========    ============   ============================================
Column      Type           Meaning
========    ============   ============================================
1           Int64/Double   OBT times
2           Double         SCET times (in ms)
3           Double         Theta angle, in radians (colatitude)
4           Double         Phi angle, in radians (longitude)
5           Double         Psi angle, in radians
========    ============   ============================================

The name of the columns is not important. The first table HDU will be
read. Since Ereline uses CFITSIO, you can use the extended name
support if the layout of the files does not fully conform to this.


Differenced data
----------------

Like pointing files, differenced data files can be provided either in
FITS or *Squeezer* format. The code assumes that each files contains
data for one OD.

Files in FITS format should have the following columns:

======= ============= ============ =============================
Column  Name          Type         Meaning
======= ============= ============ =============================
1       ``OBT``       Int64/Double OBT times
2       ``SCET``      Double       SCET times (in ms)
3       ``SKYLOAD``   Double       Sky/load difference, in Volt
4       ``FLAGS``     Int32        Flags
======= ============= ============ =============================

The column names specified in the table are used when new differenced
files are saved. They are not used during loading.


Spacecraft velocity
-------------------

The file containing the velocity of the spacecraft as a function of
time must be a FITS file whose first table HDU contains the following
columns:

========    ============   ============================================
Column      Type           Meaning
========    ============   ============================================
1           Double         SCET times (in ms)
2           Double         Velocity along the X axis, in m/s
3           Double         Velocity along the Y axis, in m/s
4           Double         Velocity along the Z axis, in m/s
========    ============   ============================================


Calibration constants
---------------------

Calibration constants are read/saved in the first HDU table of FITS
files. The following columns are used:

========    ============   ============================================
Column      Type           Meaning
========    ============   ============================================
1           Int32          Pointing ID
2           Double         Calibration constant, in Kelvin/Volt
3           Double         Offset, in Kelvin
========    ============   ============================================
