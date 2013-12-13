{
    "common": {
        "od_range": {
            "first": 91,
            "last": 563
        },
        
        // These are the radiometers we want to analyze
        "radiometer": ["LFI18M", "LFI18S", "LFI23M", "LFI23S"],

        /* This is the template to use for the file names containing
           the beams. The directory in which they are kept is
           specified in a separate JSON file, see
           "data_storage_layout" below. */
        "beam_alm": "beam_lfi{horn}_{arm}.alm",

        /* This .json file describes how standard data files (e.g.
           pointings) should be retrieved */
        "data_storage_layout": "/datastor/ereline/data/dx10/layout.json",

        /* This can be referred later in strings using
           {base_output_dir}. Note that we do not specify the user's
           name, but rather the value of the variable "user". */
        
        "base_output_dir": "/datastor/{user}/my_tests/dx10",

        "log_file": "{base_output_dir}/ereline.log"
    },

    "apply_r": {
        "run": true,
        "mask": null,
        "output_gmf": "{base_output_dir}/gmf/LFI{horn}{arm}_{od}_gmf.fits",
        "output_tod": "{base_output_dir}/datadiff/LFI{horn}{arm}_{od}_datadiff.fits"
    },

    "dipole_fit": {
        "run": true,
        "monopole_K_CMB": 2.725,
        "solar_dipole": {
            "theta_ecl_rad": 1.765248346,
            "phi_ecl_rad": 2.995840906,
            "speed_m_s": 369000.0
        },
        "total_convolve_order": 9,
        "mask": null,
        "output_gains": "{base_output_dir}/gains/LFI{horn}{arm}_{od}_gains.fits",
        "output_fit": "{base_output_dir}/gains/LFI{horn}{arm}_{od}_fit.fits"
    },

    "da_capo": {
        "run": true,
        "constraint": true,
        "output_gains": "{base_output_dir}/da_capo/LFI{horn}{arm}_{od}_gains.fits",
        "debug": {
            "dump_solution_every_n_iterations": 5,
            "dump_file_name": "{base_output_dir}/da_capo/debug/fit_LFI{horn}{arm}_{od}_step{step}.fits"
        }
    },

    "smooth_gains": {
        "run": true,
        "fast_variations": {
            "window_length": null,
            "percent": null
        },
        "smooth_window": {
            "dipole_min": 2400,
            "dipole_max": 600,
            "window_length": null
        },
        "slow_var_percentile": null,
        "dipole_range": {
            "min_value": null,
            "max_value": null
        }
    }
}

