#!/bin/bash -x
mpixlC -qsmp=omp parallel_integral_mpi.cpp parallel_integral.cpp -o par
