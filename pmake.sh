#!/bin/bash -x
xlc++_r -qsmp=omp parallel_integral_omp.cpp parallel_integral_openmp.cpp
