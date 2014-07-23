#include "dipole_fit_results.hpp"
#include "logging.hpp"
#include <assert.h>
#include <mpi.h>

void
merge_dipole_fits_from_mpi_processes(Dipole_fit_results_t & results)
{
    Logger * log = Logger::get_instance();
    const int mpi_size = MPI::COMM_WORLD.Get_size();
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();

    // Retrieve the gain table from all the other MPI processes and put
    // them together in results.gain_table
    results.gain_table.mergeResults();

    // Since odd and even MPI processes work on different radiometer arms
    // (M/S), we discard those gains that do not belong to the arm
    // being analyzed by the current MPI process
    outDipFitTable.selectRadiometerGains(mpi_rank % 2, 2, idsRange);
}
