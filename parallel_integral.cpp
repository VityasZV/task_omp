#include "parallel_integral.hpp"

#include <mpi.h>
#include <iostream>
#include <cmath>

namespace parallel_integral {

	namespace {
        	extern "C" double omp_get_wtime(void);

		static const double kAccuracy = 0.0001;


		double Function(const double& arg) {
			return 3*std::pow(arg, -1) + 2;
		}

		struct AccuracyParameters {
			unsigned long parts;
			double step;
			const Limits limits;
			AccuracyParameters(std::ifstream& input_file): limits(input_file){
				parts = 2;
				step = (limits.right - limits.left) / parts;
			}
			void Increment() {
				parts *= 2;
				step = (limits.right - limits.left) / parts;
			}
			operator double() {
				return step;
			}
		};

	} //namespace

	Limits::Limits(std::ifstream& input_file) {
		input_file.open("input.txt");
		input_file >> left >> right;
		input_file.close();
	}

	ResultAndTime ComputeIntegral(int *argc_ptr, char ***argv_ptr) {
		std::ifstream input_file;
		AccuracyParameters accuracy_parameters(input_file);
		double previous_result = 0, result = 0;
		unsigned long i = 0;
        double time = omp_get_wtime();
		double mpi_time = MPI_Wtime();

		int ierr = MPI_Init(argc_ptr, argv_ptr);//initializing mpi
		int amount_of_processes, process_id;
		if (ierr != MPI_SUCCESS) {
			return ResultAndTime(-1, -1);//костыли потому что не получилось указать компилятору использование std::exception
		}

		ierr = MPI_Comm_size(MPI_COMM_WORLD, &amount_of_processes); //get number of processes
		ierr = MPI_Comm_rank(MPI_COMM_WORLD, &process_id); //get the individual process_id

		do {
			previous_result = result;
			result = 0;
			accuracy_parameters.Increment();
#pragma omp parallel shared(accuracy_parameters) reduction(+:result)
			{
				#pragma omp for
				for (i = 0; i < accuracy_parameters.parts; ++i) {
					double x = accuracy_parameters.limits.left + i * accuracy_parameters.step;
					result += accuracy_parameters.step * Function(x + accuracy_parameters.step / 2);
				}
			}

		} while (fabs(result - previous_result) >= kAccuracy);

		MPI_Finalize(); // Finalize the MPI environment.

        time = omp_get_wtime() - time;
		mpi_time = MPI_Wtime() - mpi_time;
        return ResultAndTime(result, time);
	}

}
