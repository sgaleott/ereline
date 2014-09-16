.. Tutorial on how to run Ereline

Tutorial
========


Compiling the code
------------------

Refer to the file README.md in the source tree for information about
how to compile Ereline.

Obtaining the data
------------------

Once the Planck 2014 data are available, put here an explanation about
how to download them.

Running the pipeline
--------------------

The Ereline code tree contains a stand-alone executable,
``dx11d_pipeline``, which runs the calibration pipeline used in the
2014 data release. It is a command-line program which requires one and
only one parameter, the name of a text file containing the parameters
to be used in the execution of the pipeline. Such file is called the
*parameter file*.

The syntax of the parameter file follows the `JSON specification`_.
The following is an example:

.. code-block:: javascript
   :linenos:

   {
       "common": {
           "first_od": 91,
           "last_od": 563,
           "radiometer": ["LFI27M", "LFI28M"],
           "data_storage_layout": "./data_storage.json",
           "base_output_dir": "/home/foo/ereline/output/",
           "log_file": "{common.base_output_dir}/ereline.{mpi_rankNNNN}.log",
           "log_level": 4
       },

       "bin_data": {
           "run": true,
           "use_pencil_beam": false,
           "monopole_K_CMB": 2.72548,
           "solar_dipole": {
               "theta_ecl_rad": 1.7656131194951572,
               "phi_ecl_rad": 2.995889600573578,
               "speed_m_s": 370082.2332
           },
           "quality_flag": 6111248,
           "galactic_pickup": "/data/galactic_straylight_dx11d_LFI{horn}{arm}.fits",
           "total_convolve_order": 9,
           "mask": "/data/mask_calib_{frequency_GHz}GHz.fits.gz",
           "output_fit": "{common.base_output_dir}/binned_data/LFI{horn}{arm}/{odNNNN}/{pointing_idNNNNNN}.fits"
       },

       "dipole_fit": {
           "run": true,
           "use_mademoiselle": true,
           "output_gains": "{common.base_output_dir}/gains/LFI{horn}{arm}_{od}_gains.fits",
           "output_fit": "{common.base_output_dir}/gains/LFI{horn}{arm}_{od}_fit.fits",
           "debug": true
       },

       "da_capo": {
           "run": true,
           "constraint": true,
           "use_mademoiselle": false,
           "output_gains": "{common.base_output_dir}/da_capo/LFI{horn}{arm}_{od}_gains.fits"
       },

       "smooth_gains": {
           "run": true,
           "fast_variations": {
               "window_length": 300,
               "percent": 0.95
           },
           "hybrid_fit_first_pid": 1,
           "hybrid_fit_last_pid": 45957,
           "smooth_window": {
               "near_dipole_min": 2400,
               "near_dipole_max": 600,
               "window_length": 300
           },
           "slow_var_percentile": 0.99,
           "dipole_range": {
               "min_value": 3.4e-3,
               "max_value": 4.0e-3
           }
       }
   }

This file references another JSON file, whose name is assumed to be
``./data_storage.json``. An example is provided here:

.. code-block:: javascript

   {
       "pointings" : {
           "base_path" : "/data/pointings/LFI{horn}{arm}",
           "file_name_mask" : "LFI{horn}{arm}_pointings_{odNNNN}.sqz"
       },

       "differenced_data" : {
           "base_path" : "/data/datadiff/LFI{horn}{arm}",
           "file_name_mask" : "LFI{horn}{arm}_datadiff_{odNNNN}.sqz"
       },

       "ucds_file_path" : "/data/ucds-dx11-delta.db",
       "spacecraft_velocity_file" : "/data/spacecraft_velocity.fits"
   }


Assuming that the file is called ``lfi.json`` and is saved in the
user's home directory, to run the 2014 pipeline using this file the
user must run the command::

    mpirun -n NN /path/to/ereline/dx11d_pipeline ~/lfi.json

where ``NN`` is the number of MPI processes to use. (Depending on your
MPI library, you might have to change ``mpirun`` with some other
command.) So far, the number of MPI processes must be **even**, as the
code always calibrates the main and side arms of each radiometers in
parallel and allocates the same number of processes for both. The
program warns the user if ``NN`` is an even number.

.. _JSON specification: http://json.org/


Dissecting the example
----------------------

The file ``lfi.json`` presented in the previous section contains a
description of the way the LFI calibration pipeline should run and
extract the results. It is a hierarchical description, where the many
pipeline modules have their parameters specified in separate sections.
Each section is enclosed within curly brackets, and it contains a set
of key/value pairs separated by commas.

The following sections of this manual explain the purpose of each
module, as well as the parameters that can be specified in the
parameter file to tune the module's run. Here we present a handful of
features of such format:

1. Keys are always strings and must be enclosed within double quotes.
2. Values can be either a number (e.g., ``first_od``), a string (e.g.,
   ``data_storage_layout``), a list (whose elements are
   comma-separated and enclosed within square brackets, e.g.,
   ``radiometer``), or a set of key/values enclosed within curly
   brackets.
3. String values can reference other strings by enclosing them within
   curly brackets. If the string is enclosed within a block (e.g., in
   the example above, ``base_output_dir`` is within the ``common``
   block), the name must be preceded by names of all the parent blocks
   separated by a dot (e.g., ``base_output_dir`` must be referred as
   ``{common.base_output_dir}``). Nested references are allowed, so
   that this example is valid and expands ``ref2`` into ``this is a
   reference to bar, I say!``:

.. code-block:: javascript

   {
      "common": {
         "foo" : "bar",
         "ref1" : "this is a reference to {common.foo}",
         "ref2" : "{common.ref1}, I say!"
      }
   }

4. An additional JSON file (which follows the same rules) is required
   to specify how input data are saved to disk. (In the example, the
   file is ``./data_storage.json``.) At the moment, the only "input
   data" used by the pipeline are:

   a. Pointing information TODs;
   b. Differenced data;
   c. A copy of the Ultra-compact Database (UCDS);
   d. A FITS file containing the velocity of the Planck spacecraft as
      a function of time.
