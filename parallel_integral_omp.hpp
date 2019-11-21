#pragma once
#include <fstream>

namespace parallel_integral {

	struct Limits {
		double left;
		double right;
		Limits(std::ifstream& input_file);
	};

	// Parallel computing of integral
	double ComputeIntegral();

}// parallel_integral