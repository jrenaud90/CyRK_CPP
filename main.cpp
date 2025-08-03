/* This file and its functions are just used for crude testing. It is not considered truely apart of CyRK's C++ backend. */
#include <iostream>
#include <fstream>
#include <chrono>
#include <cstring>

#include <cmath>
#include "cysolve.hpp"

static void test_accur_diffeq(double* dy_ptr, double time, double* y_ptr, char* args_ptr, PreEvalFunc pre_eval_func)
{
    const double y0 = y_ptr[0];
    const double y1 = y_ptr[1];

    dy_ptr[0] = std::sin(time) - y1;
    dy_ptr[1] = std::cos(time) + y0;
}

static void test_diffeq(double* dy_ptr, double time, double* y_ptr, char* args_ptr, PreEvalFunc pre_eval_func)
{
    const double y0 = y_ptr[0];
    const double y1 = y_ptr[1];

    const double e1 = (1. - 0.01 * y1);
    const double e2 = (0.02 * y0 - 1.);

    dy_ptr[0] = e1 * y0;
    dy_ptr[1] = e2 * y1;
}

static void test_extra_diffeq(double* dy_ptr, double time, double* y_ptr, char* args_ptr, PreEvalFunc pre_eval_func)
{
    const double y0 = y_ptr[0];
    const double y1 = y_ptr[1];

    const double e1 = (1. - 0.01 * y1);
    const double e2 = (0.02 * y0 - 1.);

    dy_ptr[0] = e1 * y0;
    dy_ptr[1] = e2 * y1;

    dy_ptr[2] = e1;
    dy_ptr[3] = e2;
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

void test_regular(
    double t_end = 500.0,
    bool save_dense = false,
    size_t len_t_eval = 0, // 7070 == 2x; 1767 == 0.5x for tspan of (0, 500))
    unsigned int num_extra = 0,
    ODEMethod method = ODEMethod::RK45,
    bool use_resets = false
)
{
    DiffeqFuncType diffeq_func = test_diffeq;
    if (num_extra > 0)
    {
        DiffeqFuncType diffeq_func = test_extra_diffeq;
    }

    double t_start = 0.0;
    int max_i = 1000;

    std::vector<double> y0_vec = std::vector<double>(2);
    y0_vec[0] = 20.0; // Initial condition for y0
    y0_vec[1] = 20.0; // Initial condition for y0

    std::string msg = "No message";

    unsigned int num_y = 2;

    std::vector<double> rtols_vec = std::vector<double>(1);
    rtols_vec[0] = 1.0e-7; // Relative tolerance for y0 and y1
    std::vector<double> atols_vec = std::vector<double>(1);
    atols_vec[0] = 1.0e-8; // Absolute tolerance for y0 and y1

    size_t expected_size = 0;

    double running_sum = 0.0;
    size_t final_size;
    size_t sol_size;
    size_t num_interps = 0;

    std::vector<double> t_eval_vec = std::vector<double>(0);
    if (len_t_eval > 0)
    {
        t_eval_vec = linspace(t_start, t_end, len_t_eval);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();
    int k = 0;
    int k_max = 15;
    double total_runner = 0.0;
    double t_at_20;
    double y_at_20[2] = { -99.0, -99.0 };
    double* y_at_20_ptr = &y_at_20[0];
    double y_int_at_20[2];
    double* y_int_at_20_ptr = &y_int_at_20[0];
    std::unique_ptr<CySolverResult> result_uptr_regular;
    std::unique_ptr<CySolverResult> result_uptr_forresets = std::make_unique<CySolverResult>(method);
    CySolverResult* result_ptr = nullptr;

    while (k < k_max)
    {
        running_sum = 0.0;

        for (size_t i = 0; i < max_i; i++)
        {

            t1 = std::chrono::high_resolution_clock::now();
            if (use_resets)
            {
                // TODO
                result_ptr = result_uptr_forresets.get();
            }
            else
            {
                result_uptr_regular = baseline_cysolve_ivp(
                    diffeq_func,
                    t_start,
                    t_end,
                    y0_vec,
                    method,
                    expected_size,
                    num_extra,
                    std::nullopt,
                    100000,
                    2000,
                    save_dense,
                    t_eval_vec,
                    std::nullopt, // pre_eval_func
                    rtols_vec,
                    atols_vec,
                    INF, // max_step_size
                    0.0 // first_step_size
                );
                result_ptr = result_uptr_regular.get();
            }

            final_size = result_ptr->size;
            msg = result_ptr->message;
            num_interps = result_ptr->num_interpolates;
            t_at_20 = result_ptr->time_domain_vec[4];
            y_at_20_ptr = &result_ptr->solution[(num_y + num_extra) * 4];
            result_ptr->call(t_at_20, y_int_at_20_ptr);

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
        std::cout << "Message: " << msg << std::endl;
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
    datastream.precision(32);

    if (result_ptr)
    {
        for (size_t i = 0; i < final_size; i++)
        {
            datastream << result_ptr->time_domain_vec[i] << ", " << result_ptr->solution[2 * i] << ", " << result_ptr->solution[2 * i + 1] << std::endl;
        }
    }
    datastream.close();
    std::cout << "Data saved." << std::endl;

    //std::cin.get();
}

int main(){
    /* This file and its functions are just used for crude testing. It is not considered truely apart of CyRK's C++ backend. */

    /* Speed Runs (JPR Work Desktop):
    * No Dense; No t_eval
    *   755.773us; 610.87us
    * 
    * With Dense; No t_eval
    *   1089.25us; 1094.09
    * 
    * With Dense; With 2x t_eval
    *   1246.61us; 1239.9
    * 
    * With Dense; With 0.5x t_eval
    *   1135.45us; 1136.49
    *
    * No Dense; With 0.5x t_eval
    *   822.59us; 820.39
    * 
    * No Dense; With 2x t_eval
    *   902.975us; 901.372
    * 
    * No Dense; With 2x t_eval; With extra on
    *   1000.79us; 952.867; 941.13
    * 
    * No Dense; With 0.5x t_eval; With extra on
    *   843.989us; 837.096
    */

    test_regular(
        500.0, // t_end
        false, // Dense
        0,     // len t_eval 7070 == 2x; 1767 == 0.5x for tspan of (0, 500))
        0,     // num extra
        ODEMethod::RK45,      // Method
        false  // use_resets
    );

    // 1333.07


    return 0;
}
