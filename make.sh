#!/bin/bash -x
mpixlC -qsmp=omp parallel_integral_mpi.cpp parallel_integral.cpp mpi_info.cpp -o par
