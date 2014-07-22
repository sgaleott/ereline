#include "configuration.hpp"
#include "logging.hpp"
#include "dipole_fit.hpp"
#include "da_capo.hpp"
#include "datatypes.hpp"
#include "smooth_gains.hpp"
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
read_ahf_info(SQLite3Connection & ucds,
	      const Configuration & program_config,
	      std::vector<Pointing_t> & list_of_pointings)
{
    int first_od, last_od;
    first_od = program_config.get<int>("common.first_od");
    last_od = program_config.get<int>("common.last_od");

    loadPointingInformation(ucds, first_od, last_od, list_of_pointings);

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

int
inner_main(int argc, const char ** argv)
{
    if(MPI::COMM_WORLD.Get_size() % 2 != 0) {
	if(MPI::COMM_WORLD.Get_rank() == 0) {
	    std::cerr << "Error: the program expects to be run on "
		"an even number of MPI processes.\n";
	}
	return 1;
    }

    if(argc != 2) {
	if(MPI::COMM_WORLD.Get_rank() == 0)
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
    SQLite3Connection ucds(ucds_file_path.c_str());

    // Convert pointing information into appropriate C++ structures
    std::vector<Pointing_t> list_of_pointings;
    read_ahf_info(ucds, program_config, list_of_pointings);

    const auto radiometer_nodes =
	program_config.ptree.get_child("common.radiometer");
    for(const auto & cur_node : radiometer_nodes) {
	LfiRadiometer radiometer(cur_node.second.get<std::string>(""));

	log->info(boost::format("Processing radiometer %1%")
		  % radiometer.shortName());

	run_dipole_fit(ucds,
		       radiometer, 
		       program_config, 
		       storage_config, 
		       list_of_pointings);
	run_da_capo(program_config, storage_config);
	run_smooth_gains(program_config, storage_config);
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
