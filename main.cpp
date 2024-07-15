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


std::vector<double> linspace(double start, double end, size_t num_in)
{

    std::vector<double> linspaced;
    linspaced.reserve(num_in);

    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for (int i = 0; i < num - 1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); // I want to ensure that start and end
    // are exactly the same as the input
    return linspaced;
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
    
    size_t expected_size = 0;
    
    double running_sum = 0.0;
    size_t final_size;
    size_t sol_size;
    size_t num_interps = 0;


    const bool dense_output = true;


    size_t len_t_eval = 7070;  // 7070 == 2x; 1767 == 0.5x for tspan of (0, 500)
    std::vector<double> t_eval = linspace(time_span_ptr[0], time_span_ptr[1], len_t_eval);
    const double* t_eval_ptr = &t_eval[0];

    //len_t_eval = 0;
    //t_eval_ptr = nullptr;

    std::shared_ptr<CySolverResult> result;

    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();
    int k = 0;
    int k_max = 20;
    double total_runner = 0.0;
    double t_at_20;
    double y_at_20[2] = { -99.0, -99.0 };
    double* y_at_20_ptr = &y_at_20[0];
    double y_int_at_20[2];
    double* y_int_at_20_ptr = &y_int_at_20[0];
    while (k < k_max)
    {
        running_sum = 0.0;

        for (size_t i = 0; i < max_i; i++)
        {

            t1 = std::chrono::high_resolution_clock::now();
            result = baseline_cysolve_ivp(
                diffeq_func,
                time_span_ptr,
                y0_ptr,
                num_y,
                method,
                expected_size,
                num_extra,
                nullptr,
                dense_output,
                t_eval_ptr,
                len_t_eval,
                100000,
                2000,
                rtol,
                atol
            );

            final_size = result->size;
            msg_ptr = result->message_ptr;
            num_interps = result->num_interpolates;
            t_at_20 = result->time_domain[4];
            y_at_20_ptr = &result->solution[2 * 4];
            result->call(t_at_20, y_int_at_20_ptr);

            //std::cout << result->message_ptr << std::endl;
            sol_size = final_size * 2;

            t2 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::micro> ms_double = t2 - t1;
            running_sum += ms_double.count();
            //std::cout << ms_double.count() << "ms\n";
        }

        std::cout << "SIZE: " << final_size << std::endl;
        std::cout << "NUM INTERPS: " << num_interps << std::endl;
        std::cout << "20 t = " << t_at_20 << "; y0 = " << y_at_20_ptr[0] << ", y1 = " << y_at_20_ptr[1] << std::endl;
        std::cout << "DENSE INTERP:: " << "y0 = " << y_int_at_20_ptr[0] << ", y1 = " << y_int_at_20_ptr[1] << std::endl;
        std::cout << "Message: " << msg_ptr << std::endl;
        std::cout << "AVERAGE: " << running_sum / max_i << "us\n" << std::endl;
        if (k > 3)
        {
            total_runner += running_sum / max_i;
        }
        k += 1;
        //break;
    }

    std::cout << "Done! Final Avg: " << total_runner / (k_max - 4) << std::endl << "Saving data..." << std::endl;

    
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
