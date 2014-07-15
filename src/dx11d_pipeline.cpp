#include "configuration.hpp"
#include "logging.hpp"
#include "dipole_fit.hpp"
#include "da_capo.hpp"
#include "smooth_gains.hpp"
#include "sqlite3xx.hpp"
#include "ahf_info.hpp"
#include "io.hpp"

#include <iostream>

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
read_configuration_and_set_up(const std::string & configuration_file_name,
			      Configuration & program_config,
			      Configuration & storage_config)
{
    program_config.read_from_json(configuration_file_name);
    storage_config.read_from_json(program_config.getWithSubst("common.data_storage_layout"));

    // Configure the logger according to the settings found in the
    // configuration file
    Logger * log = Logger::get_instance();
    program_config.configure_logging();
}

////////////////////////////////////////////////////////////////////////////////

void
read_ahf_info(SQLite3Connection & ucds,
	      const Configuration & program_config,
	      std::vector<Pointing_t> & list_of_pointings,
	      std::vector<Od_t> & list_of_ods)
{
    loadPointingInformation(ucds,
			    program_config.get<int>("common.first_od"),
			    program_config.get<int>("common.last_od"),
			    list_of_pointings);
    list_of_ods = build_od_list(list_of_pointings);
}

////////////////////////////////////////////////////////////////////////////////

int
main(int argc, const char ** argv)
{
    if(argc != 2) {
	print_help();
	return 1;
    }

    Configuration program_config;
    Configuration storage_config;
    read_configuration_and_set_up(argv[1], program_config, storage_config);

    Logger * log = Logger::get_instance();

    std::string ucds_file_path = storage_config.getWithSubst("ucds_file_path");
    log->info(boost::format("UCDS file path: %1%") % ucds_file_path);
    SQLite3Connection ucds(ucds_file_path.c_str());

    std::vector<Pointing_t> list_of_pointings;
    std::vector<Od_t> list_of_ods;

    read_ahf_info(ucds, program_config, list_of_pointings, list_of_ods);
    log->info(boost::format("Number of pointings to process: %1% (%2%-%3%), number of "
			    "ODs to process: %4% (%5%-%6%)")
	      % list_of_pointings.size()
	      % list_of_pointings.front().id
	      % list_of_pointings.back().id
	      % list_of_ods.size()
	      % list_of_ods.front().od
	      % list_of_ods.back().od);

    run_dipole_fit(program_config);
    run_da_capo(program_config);
    run_smooth_gains(program_config);

    return 0;
}
