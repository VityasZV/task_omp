#include "mpi_info.hpp"

#include <mpi.h>
#include <exception>

namespace mpi_info {
    MPI::MPI(int *argc_ptr, char ***argv_ptr){
        int ierr = MPI_Init (argc_ptr, argv_ptr);
        //  Initialize MPI.
        if (ierr != 0)
        {
            throw std::runtime_error("MPI_Init returned nonzero ierr.");
        }
    }
    MPI::~MPI(){
        MPI_Finalize();
    }
}// mpi_info

