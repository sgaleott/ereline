#include "configuration.hpp"
#include "logging.hpp"
#include "dipole_fit.hpp"
#include "dipole_fit_results.hpp"
#include "da_capo.hpp"
#include "da_capo_results.hpp"
#include "datatypes.hpp"
#include "mademoiselle.hpp"
#include "smooth_gains.hpp"
#include "smooth_gains_results.hpp"
#include "sqlite3xx.hpp"
#include "ahf_info.hpp"
#include "io.hpp"

#include <iostream>

#include <mpi.h>

////////////////////////////////////////////////////////////////////////////////

void
print_help()
{
    std::cout << "Usage:\n"
              << "\n"
              << "    dx11d_pipeline CONFIGURATION_FILE\n"
              << "\n"
              << "Run the DX11\u03b4 pipeline using the parameters found in\n"
              << "CONFIGURATION_FILE (a JSON file).\n"
              << "\n"
              << "For further information, refer to the site\n"
              << "\n"
              << "    http://belen:8080/tomasi/ereline\n"
              << "\n"
              << "(password protected).\n";
}

////////////////////////////////////////////////////////////////////////////////

void
set_mpi_variables_in_configuration(int mpi_rank,
                                   int mpi_size,
                                   Configuration & conf)
{
    conf.set_variable("mpi_rank",
                      (boost::format("%d") % (mpi_rank + 1)).str());
    conf.set_variable("mpi_rankNN",
                      (boost::format("%02d") % (mpi_rank + 1)).str());
    conf.set_variable("mpi_rankNNNN",
                      (boost::format("%04d") % (mpi_rank + 1)).str());

    conf.set_variable("mpi_size",
                      (boost::format("%d") % mpi_size).str());
    conf.set_variable("mpi_sizeNN",
                      (boost::format("%02d") % mpi_size).str());
    conf.set_variable("mpi_sizeNNNN",
                      (boost::format("%04d") % mpi_size).str());
}

////////////////////////////////////////////////////////////////////////////////

void
read_configuration_and_set_up(const std::string & configuration_file_name,
                              Configuration & program_config,
                              Configuration & storage_config)
{
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
    const int mpi_size = MPI::COMM_WORLD.Get_size();

    set_mpi_variables_in_configuration(mpi_rank, mpi_size, program_config);
    set_mpi_variables_in_configuration(mpi_rank, mpi_size, storage_config);

    program_config.read_from_json(configuration_file_name);
    storage_config.read_from_json(program_config.getWithSubst("common.data_storage_layout"));

    // Configure the logger according to the settings found in the
    // configuration file
    program_config.configure_logging();

    Logger * log = Logger::get_instance();
    log->set_mpi_rank(mpi_rank, mpi_size);
}

////////////////////////////////////////////////////////////////////////////////

void
read_ahf_info(Sqlite_connection_t & ucds,
              const Configuration & program_config,
              std::vector<Pointing_t> & list_of_pointings)
{
    int first_od, last_od;
    first_od = program_config.get<int>("common.first_od");
    last_od = program_config.get<int>("common.last_od");

    load_pointing_information(ucds, first_od, last_od, list_of_pointings);

    Logger * log = Logger::get_instance();
    log->info(boost::format("Number of pointings to process: %1% (%2%-%3%), "
                            "ODs to process are in the range %4% - %5%")
              % list_of_pointings.size()
              % list_of_pointings.front().id
              % list_of_pointings.back().id
              % first_od
              % last_od);
}

////////////////////////////////////////////////////////////////////////////////

const std::string
dipole_fit_name(const Configuration & program_config,
        const Lfi_radiometer_t & radiometer,
        bool da_capo)
{
    if(da_capo) {
    return (boost::format("%s/da_capo_dipole_fit_table_%s.fits")
        % program_config.getWithSubst("common.base_output_dir")
        % radiometer.shortName()).str();
    } else {
    return (boost::format("%s/dipole_fit/dipole_fit_table_%s.fits")
        % program_config.getWithSubst("common.base_output_dir")
        % radiometer.shortName()).str();
    }
}

////////////////////////////////////////////////////////////////////////////////

int
inner_main(int argc, const char ** argv)
{
    const int mpi_size = MPI::COMM_WORLD.Get_size();
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();

    if(mpi_size % 2 != 0) {
        if(mpi_rank == 0) {
            std::cerr << "Error: the program expects to be run on "
                "an even number of MPI processes.\n";
        }
        return 1;
    }

    if(argc != 2) {
        if(mpi_rank == 0)
            print_help();

        return 1;
    }

    // Load the JSON files containing the configuration and set up a
    // few structures
    Configuration program_config;
    Configuration storage_config;
    read_configuration_and_set_up(argv[1], program_config, storage_config);

    Logger * log = Logger::get_instance();

    // Open a connection with the UCDS database
    std::string ucds_file_path = storage_config.getWithSubst("ucds_file_path");
    log->info(boost::format("UCDS file path: %1%") % ucds_file_path);
    Sqlite_connection_t ucds(ucds_file_path.c_str());

    // Convert pointing information into appropriate C++ structures
    std::vector<Pointing_t> list_of_pointings;
    read_ahf_info(ucds, program_config, list_of_pointings);

    const auto radiometer_nodes =
        program_config.ptree.get_child("common.radiometer");
    for(const auto & cur_node : radiometer_nodes) {
        Lfi_radiometer_t radiometer(cur_node.second.get<std::string>(""));

        log->info(boost::format("Processing radiometer %1%")
                  % radiometer.shortName());

        Dipole_fit_results_t dipole_fit_results;
        run_dipole_fit(ucds,
                       radiometer,
                       program_config,
                       storage_config,
                       list_of_pointings,
                       dipole_fit_results);
        if(mpi_rank < 2) {
            save_dipole_fit_results(ensure_path_exists(dipole_fit_name(program_config,
                                                                       radiometer, false)),
                                    radiometer, dipole_fit_results);
        }

        Da_capo_results_t da_capo_results;
        if(program_config.get<bool>("da_capo.use_mademoiselle")) {
            run_mademoiselle(program_config,
                             storage_config,
                             radiometer,
                             list_of_pointings,
                             dipole_fit_results,
                             da_capo_results);
        } else {
            run_da_capo(program_config,
                        storage_config,
                        radiometer,
                        list_of_pointings,
                        dipole_fit_results,
                        da_capo_results);
        }

        if(mpi_rank < 2) {
            save_dipole_fit_results(ensure_path_exists(dipole_fit_name(program_config,
                                                                       radiometer, true)),
                                    radiometer, dipole_fit_results);
        }

        Smooth_gains_results_t smooth_results;
        run_smooth_gains(ucds,
                         program_config,
                         storage_config,
                         radiometer,
                         list_of_pointings,
                         dipole_fit_results,
                         da_capo_results,
                         smooth_results);
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int
main(int argc, const char ** argv)
{
    MPI::Init();
    int return_code = inner_main(argc, argv);
    MPI::Finalize();

    return return_code;
}
