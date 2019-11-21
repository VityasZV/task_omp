// parallel_integral_openmp.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <omp.h>
#include "parallel_integral_omp.hpp"

int main()
{
	std::cout << parallel_integral::ComputeIntegral();
}
