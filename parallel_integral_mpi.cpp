// parallel_integral_openmp.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <ctime>
#include "parallel_integral.hpp"

int main(int argc, char **argv)
{
    parallel_integral::ResultAndTime result_and_time = parallel_integral::ComputeIntegral(&argc, &argv);
    std::cout << "result = " << result_and_time.result;
    std::cout << "\ntime = " << result_and_time.time << std::endl;
    return 0;
}
