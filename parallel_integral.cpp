#include "parallel_integral.hpp"

#include "mpi_info.hpp"
#include <mpi.h>

#include <iostream>
#include <cmath>
#include <stdio.h>
#include <cstdlib>

namespace parallel_integral {

	namespace {
        extern "C" double omp_get_wtime(void); //possibly unneeded in MPI

		static const double kAccuracy = 0.000001;


		double Function(const double& arg) {
			return 1;
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
		double mpi_time = MPI_Wtime();
		double *collecting_result = NULL; //buffer for master process
		MPI_Request* requests = NULL;
		MPI_Status* statuses = NULL;
		mpi_info::MPI mpi_statistics(argc_ptr, argv_ptr);
		if (mpi_statistics.ierr != MPI_SUCCESS) {
			return ResultAndTime(-1, -1); //костыли потому что не получилось указать компилятору использование std::exception
		}
		//
        //  Get the number of processes.
        //
        mpi_statistics.ierr = MPI_Comm_size(MPI_COMM_WORLD, &mpi_statistics.amount_of_processes);
        //
        //  Get the individual process ID.
        //
        mpi_statistics.ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_statistics.process_id);
		/*
			Giving pieces from accuracy_parameters.parts for each process,
			and then in for-cycle each process is using threads for computation of sum. 
		*/
		if (mpi_statistics.process_id == 0){
			std::cout << "PROCESS: "<<mpi_statistics.process_id <<" - allocates memory for arrays" << std::endl;
			requests = (MPI_Request*)std::malloc(sizeof(MPI_Request) * mpi_statistics.amount_of_processes);
			statuses = (MPI_Status*)std::malloc(sizeof(MPI_Status) * mpi_statistics.amount_of_processes);
			collecting_result = (double*)std::malloc(sizeof(double) * mpi_statistics.amount_of_processes); //making an array for results coming from each process
		}
		do {
			if (mpi_statistics.process_id == 0){
				std::cout << "PROCESS: "<<mpi_statistics.process_id << "changes params and broadcast" << std::endl;
				previous_result = result;
				result = 0;
				accuracy_parameters.Increment();
				double temp = 0
				//cycle that is instead of broadcast
				for (int k = 1; k < mpi_statistics.amount_of_processes){
					MPI_Send(&temp, 1, MPI_DOUBLE, k, 0, MPI_COMM_WORLD);
				}
				//commented next line bacause it is possibly wrong
				//MPI_Bcast(&result, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // sending message for the rest processes, all must wait
			}
			else {
				double temp = 0;
				MPI_Recv(&temp, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				std::cout << "PROCESS: "<<mpi_statistics.process_id << "received broadcast";
				//now I'm sure that every process waited for master 
			}
			double result_in_process = 0;
			unsigned long for_work = accuracy_parameters.parts / mpi_statistics.amount_of_processes;
			unsigned long remains = accuracy_parameters.parts % mpi_statistics.amount_of_processes;
#pragma omp parallel shared(accuracy_parameters) reduction(+:result_in_process)
			{
				unsigned long start = mpi_statistics.process_id * for_work;
				unsigned long finish = start + for_work + mpi_statistics.process_id == 0 ? remains : 0;
				#pragma omp for
				for (i = start; i < finish; ++i) {
					double x = accuracy_parameters.limits.left + i * accuracy_parameters.step;
					result_in_process += accuracy_parameters.step * Function(x + accuracy_parameters.step / 2);
				}

			}
			collecting_result[mpi_statistics.process_id] = result_in_process;
			if (mpi_statistics.process_id != 0){
				std::cout << "PROCESS: "<<mpi_statistics.process_id << "sending nothing to master" << std::endl;
				MPI_Send(NULL, 0, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); //sending nothing but master should wait for others
			}
			if (mpi_statistics.process_id == 0) {
				int count = 0;
				std::cout << "PROCESS: "<<mpi_statistics.process_id << "is waiting for others finally" << std::endl; //here was infinite loop last time

				MPI_Recv(NULL, 0, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				count+=1;

				std::cout << "PROCESS: "<<mpi_statistics.process_id << "waited for count=" << count << std::endl;
				//MPI_Waitall(mpi_statistics.amount_of_processes - 1, &requests[1], &statuses[1]); // дождался всех 
				std::cout << "DEBUG: Master waited for all " << mpi_statistics.amount_of_processes - 1 << std::endl;
				for (int j = 0; j < mpi_statistics.amount_of_processes; ++j){
					result += collecting_result[j]; // master is collecting data from others and forming a result
				}
			}

		} while (fabs(result - previous_result) >= kAccuracy);
        time = omp_get_wtime() - time;
		mpi_time = MPI_Wtime() - mpi_time;
		std::free(requests);
		std::free(statuses);
		std::free(collecting_result);
        return ResultAndTime(result, time);
	}

}
