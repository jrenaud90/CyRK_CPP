#include <iostream>
#include <fstream>
#include <chrono>
#include <cstring>

#include "cysolve.hpp"

static void test_diffeq(double* dy_ptr, double time, double* y_ptr, const double* args_ptr)
{
    const double y0 = y_ptr[0];
    const double y1 = y_ptr[1];

    dy_ptr[0] = (1. - 0.01 * y1) * y0;
    dy_ptr[1] = (0.02 * y0 - 1.) * y1;
}

int main(){
    //std::cout << "Running..." << std::endl;

    DiffeqFuncType diffeq_func = test_diffeq;

    double time_span[2] = {0.0, 500.0};
    double* time_span_ptr = &time_span[0];
    unsigned int method = 1;
    int max_i = 1000;
    
    double y0[2] = {20.0, 20.0};
    double* y0_ptr = &y0[0];
    char msg[512];
    char* msg_ptr = &msg[0];

    unsigned int num_y = 2;
    unsigned int num_extra = 0;

    double rtol = 1.0e-7;
    double atol = 1.0e-8;
    
    size_t expected_size = 300;
    
    double running_sum = 0.0;
    size_t final_size;
    size_t sol_size;
    std::shared_ptr<CySolverResult> result;
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();
    int k = 0;
    double total_runner = 0.0;
    while (k < 20)
    {
        running_sum = 0.0;

        for (size_t i = 0; i < max_i; i++)
        {

            t1 = std::chrono::high_resolution_clock::now();
            result = cysolve_ivp(
                diffeq_func,
                time_span_ptr,
                y0_ptr,
                num_y,
                method,
                expected_size,
                num_extra,
                nullptr,
                100000,
                2000,
                rtol,
                atol
            );

            final_size = result->size;
            //std::cout << result->message_ptr << std::endl;
            sol_size = final_size * 2;

            t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::micro> ms_double = t2 - t1;
            running_sum += ms_double.count();
            //std::cout << ms_double.count() << "ms\n";
        }

        std::cout << "SIZE: " << final_size << std::endl;
        std::cout << "Message: " << msg_ptr << std::endl;
        std::cout << "AVERAGE: " << running_sum / max_i << "us\n" << std::endl;
        total_runner += running_sum / max_i;
        k += 1;
        //break;
    }

    std::cout << "Done! Final Avg: " << total_runner / 20 << std::endl << "Saving data..." << std::endl;

    
    std::ofstream datastream;

    datastream.open("out.dat");
    if (!datastream)
    {
        std::cout << "Could not open data file." << std::endl;
    }
    datastream.precision(6);

    for (size_t i = 0; i < final_size; i++)
    {
        datastream << result->time_domain[i] << ", " << result->solution[2 * i] << ", " << result->solution[2 * i + 1] << std::endl;
    }
    datastream.close();
    std::cout << "Data saved." << std::endl;

    //std::cin.get();

    return 0;
}