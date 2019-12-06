#include "parallel_integral.hpp"

#include "mpi_info.hpp"

#include <iostream>
#include <cmath>

namespace parallel_integral {

	namespace {
        extern "C" double omp_get_wtime(void); //possibly unneeded in MPI

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
        double time = omp_get_wtime(); //needs change to MPI
		mpi_info::MPI mpi_statistics(argc_ptr, argv_ptr);
		int ierr;
		//
        //  Get the number of processes.
        //
        ierr = MPI_Comm_size(MPI_COMM_WORLD, &mpi_statistics.amount_of_processes);
        //
        //  Get the individual process ID.
        //
        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_statistics.process_id);
        //
        //  Process 0 prints an introductory message.
        //
		do {
			previous_result = result;
			result = 0;
			accuracy_parameters.Increment();
			{
				for (i = 0; i < accuracy_parameters.parts; ++i) {
					double x = accuracy_parameters.limits.left + i * accuracy_parameters.step;
					result += accuracy_parameters.step * Function(x + accuracy_parameters.step / 2);
				}
			}

		} while (fabs(result - previous_result) >= kAccuracy);
        time = omp_get_wtime() - time;
        return ResultAndTime(result, time);
	}

}
