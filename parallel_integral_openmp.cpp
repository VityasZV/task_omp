// parallel_integral_openmp.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <omp.h>
#include <ctime>
#include "parallel_integral_omp.hpp"

int main()
{
    parallel_integral::ResultAndTime result_and_time = parallel_integral::ComputeIntegral();
    std::cout << "result = " << result_and_time.result;
    std::cout << "\ntime = " << result_and_time.time;
}
