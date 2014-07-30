#include <mpi.h>

#include "dipole_fit_results.hpp"
#include "fits_object.hpp"
#include "io.hpp"
#include "logging.hpp"

////////////////////////////////////////////////////////////////////////////////

void
save_dipole_fit_results(const std::string & file_name,
                        const Lfi_radiometer_t & radiometer,
                        const Dipole_fit_results_t & results,
                        const std::string & comment)
{
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
    const int mpi_size = MPI::COMM_WORLD.Get_size();

    FitsObject file;
    file.create(file_name, true);

    // 1. Gain table
    save_gain_table(file, radiometer, results.gain_table, comment);

    // 2. Fits
    {
        std::vector<fitscolumn> columns {
            { "PID", "", fitsTypeC<int>(), 1 },
            { "GAINV", "", fitsTypeC<double>(), 1 },
            { "OFFSET", "", fitsTypeC<double>(), 1 },
            { "DIPMIN", "", fitsTypeC<double>(), 1 },
            { "DIPMAX", "", fitsTypeC<double>(), 1 }};
        file.insertBINtable(columns, "DIPFITS");

        for(size_t idx = 0; idx < results.list_of_fits.size(); ++idx) {
            const Dipole_fit_t &fit(results.list_of_fits[idx]);
            const size_t fits_row = idx + 1;
            file.writeElement(1, fit.pointingID, fits_row);
            file.writeElement(2, fit.gainv, fits_row);
            file.writeElement(3, fit.offset, fits_row);
            file.writeElement(4, fit.minDipole, fits_row);
            file.writeElement(5, fit.maxDipole, fits_row);
        }
    }

    // 3. Mask
    {
        std::vector<fitscolumn> columns {
            { "VALUE", "", fitsTypeC<float>(), 1 } };
        file.insertBINtable(columns, "MASK");

        file.writeColumn(1, results.mask.pixels);
        file.setKey("NSIDE", results.mask.nside);
        file.setKey("ORDERING", std::string("NEST"));
    }

    // 4. PIDs per process
    {
        std::vector<fitscolumn> columns {
            { "NUM_PIDS", "", fitsTypeC<int>(), 1 } };
        file.insertBINtable(columns, "MPI_PIDS");

        file.writeColumn(1, results.pids_per_process);
        file.setKey("MPI_SIZE", mpi_size);
    }
}
