// parallel_integral_openmp.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <omp.h>
#include "parallel_integral_omp.hpp"

int main()
{
    result_and_time = parallel_integral::ComputeIntegral();
    std::cout << "result = " << result_and_time.result;
    std::cout << "time = " << result_and_time.time;
}
