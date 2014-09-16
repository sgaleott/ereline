.. Overview of the modules used by Ereline

Overview of the modules used by Ereline
=======================================

This section describes the way the modules in the pipeline are used by
the ``dx11d_pipeline`` program to produce a set of calibration
constants from a set of inputs (pointings, differenced data,
spacecraft velocities, etc.). The explanation shows how each module is
called by ``dx11d_pipeline`` (a C++11 program), so that interested
users can develop their own pipeline by combining the modules in
different ways.


Projecting TODS into rings
--------------------------

(In the nominal version of the pipeline developed by the LFI DPC, this
module was integrated into the dipole fitting module, which is here
described in the following section. It has been separated because in
this way the binned rings saved by this module can be reused in
subsequent runs, thus speeding up the code.)

The first step taken by the ``dx11d_pipeline`` program is to reduce
the amount of data to use for the calibration by projecting each
pointing period on a ring in a Healpix map.

The module can be run by calling the following function:

.. cpp:function:: void run_data_binning(Sqlite_connection_t & ucds, \
                                        const Lfi_radiometer_t & rad, \ 
                                        Configuration & program_conf, \
                                        Configuration & storage_conf, \
                                        const std::vector<Pointing_t> & list_of_pointings, \
                                        Data_binning_results_t & result)

   :param ucds: A SQLite3 connection to the UCDS.
   :param rad: The radiometer to process.
   :param program_conf: A :cpp:class:`Configuration` variable pointing
                        to the JSON parameter file
   :param storage_conf: A :cpp:class:`Configuration` variable pointing
                        to the JSON file describing where input data
                        files are stored.
   :param list_of_pointings: List of the pointing IDs to process
   :param result: The variable that will contain the results of the
                  computation.

   Project the samples of TODs that fall within each pointing period
   into a Healpix map. At the end, the maps (one per pointing period)
   will be saved in a compressed format. Each map contain the binned
   data as well as the simulated signal to be used by the dipole
   fitting routines (see below).

   A number of parameters used for the binning process are taken from
   ``program_conf``: this argument is not ``const``, as the ``od``
   variable needs to be updated accordingly during the execution of
   ``run_data_binning``.

   The function will make use of all the MPI processes available. If
   the execution succeeds, then each MPI process will have the full
   range of results in the ``result`` variable for one of the two
   radiometers (either M or S): for processes with even rank, the data
   in ``result`` will refer to the radiometer specified by
   ``radiometer``, for processes with odd rank ``result`` will refer
   to the twin radiometer.


Parameters in the configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+--------------------------+------------+---------------------------------------+
| Parameter                | Example    | Description                           |
+==========================+============+=======================================+
| ``run``                  | ``true``   | If ``false``, the module will not be  |
|                          |            | run. Instead, results will be loaded  |
|                          |            | from disk.                            |
+--------------------------+------------+---------------------------------------+
| ``use_pencil_beam``      | ``false``  | If ``true``, a pencil beam will be    |
|                          |            | used to produce the simulated signal. |
|                          |            | (The default is to use the 4π beam to |
|                          |            | model the dipole.)                    |
+--------------------------+------------+---------------------------------------+
| ``monopole_K_CMB``       | 2.72548    | Thermodynamic temperature of the CMB  |
|                          |            | monopole, used to produce the         |
|                          |            | simulated signal.                     |
+--------------------------+------------+---------------------------------------+
| ``solar_dipole``         | 2.72548    | A block with the keys                 |
|                          |            | ``theta_ecl_rad``, ``phi_ecl_rad``    |
|                          |            | (coordinates of the dipole axis, in   |
|                          |            | radians), and ``speed_m_s`` (speed of |
|                          |            | the Sun with respect to the CMB, in   |
|                          |            | m/s).                                 |
+--------------------------+------------+---------------------------------------+
| ``quality_flag``         | 6111248    | Integer number representing the mask  |
|                          |            | of bits in the ``FLAGS`` column of the|
|                          |            | TODs that must **not** be set for the |
|                          |            | sample to be used in the binning.     |
+--------------------------+------------+---------------------------------------+
| ``galactic_pickup``      | "gsl.fits" | Path to a FITS file containing the    |
|                          |            | ringsets for the convolution between  |
|                          |            | the beam and the Galactic signal, used|
|                          |            | to model the simulated signal. (This  |
|                          |            | parameter is used even if the flag    |
|                          |            | ``use_pencil_beam`` is ``true``.)     |
+--------------------------+------------+---------------------------------------+
| ``total_convolve_order`` | 9          | Interpolation order used for the      |
|                          |            | ringset ``galactic_pickup``.          |
+--------------------------+------------+---------------------------------------+
| ``mask``                 | "mask.fits"| Path to a Healpix map used to mask    |
|                          |            | the pixels on the Galactic plane      |
|                          |            | (zero means that the pixel will       |
|                          |            | **not** be used).                     |
+--------------------------+------------+---------------------------------------+
| ``output_fit``           |            | Where to save each compressed map.    |
|                          |            | The variable ``{pointing_id}`` should |
|                          |            | be used, either in the file name or   |
|                          |            | in the path, in order to prevent each |
|                          |            | save from overwriting the previous    |
|                          |            | one.                                  |
+--------------------------+------------+---------------------------------------+


Results
~~~~~~~

.. cpp:class:: Data_binning_results_t

   This class contains the result of a run of
   :cpp:func:`run_data_binning`. The public members of this class are
   listed in the following table:

================= ========================================== ===============================================================
Name              Type                                       Description
================= ========================================== ===============================================================
radiometer        :cpp:type:`Lfi_radiometer_t`               The radiometer which acquired the data used in the computation
binned_pids       ``std::vector<Binned_data_t>``             The compressed Healpix maps
mask              ``Healpix::Map_t<float>``                  The mask used in the computation
mpi_size          int                                        Number of MPI processes used in the computation
pids_per_process  ``std::vector<int>``                       Number of pointings processed by each MPI process
================= ========================================== ===============================================================

Example
~~~~~~~

The following code shows how to use :cpp:func:`run_data_binning` to
bin the data for radiometer LFI27M and LFI27S in the OD range 91−563:

.. code-block:: c++

   const Lfi_radiometer_t radiometer("LFI27M");

   // Open a connection to the UCDS
   Sqlite_connection_t ucds("/path/to/ucds.db");

   Configuration program_config;
   Configuration storage_config;

   // Set up program_config and storage_config
   // ...

   // Process all the data
   std::vector<Pointing_t> list_of_pointings;
   load_pointing_information(ucds, 91, 563, list_of_pointings);

   Data_binning_results_t binned_data;
   run_data_binning(ucds,
                    radiometer,
                    program_config,
                    storage_config,
                    list_of_pointings,
                    binned_data);

   // Save results for LFI27M (rank == 0) and for LFI27S (rank == 1)
   if(MPI::COMM_WORLD.Get_rank() == 0 ||
      MPI::COMM_WORLD.Get_rank() == 1)
       binned_data.save_to_disk(program_config);


Fitting the projected TODs with a model dipole
----------------------------------------------

Once the TODs have been projected into Healpix maps, it is possible to
run the dipole fitting routine :cpp:func:`run_dipole_fit`. This
function obtains a rough estimate of the calibration constant (in K/V)
by fitting each map against a model of the sky signal which includes
the dipole and the sidelobe pickup of the Galactic signal. (The
Galactic signal is assumed to have been already masked out in the
projection, see above.)


Running Da Capo
---------------


Running the smoother
--------------------
