#include <iostream>

#include "cysolver.hpp"

#include <fstream>

#include <chrono>

static void test_diffeq(double* dy_ptr, double time, double* y_ptr, double* args_ptr)
{
    const double y0 = y_ptr[0];
    const double y1 = y_ptr[1];

    dy_ptr[0] = (1. - 0.01 * y1) * y0;
    dy_ptr[1] = (0.02 * y0 - 1.) * y1;
}

int main()
{
    //std::cout << "Running..." << std::endl;

    DiffeqFuncType diffeq_func = test_diffeq;

    double time_span[2] = {0.0, 500.0};
    double* time_span_ptr = &time_span[0];
    
    double t_start = time_span[0];
    double t_end = time_span[1];

    double y0[2] = {20.0, 20.0};
    double* y0_ptr = &y0[0];

    size_t num_y = 2;

    double rtol = 1.0e-7;
    double atol = 1.0e-8;
    
    size_t expected_size = 300;
    double running_sum = 0.0;
    size_t final_size;
    size_t sol_size;
    CySolverResult* result;
    RK45* solver;
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();
    bool step_success, success;

    while (true)
    {
        running_sum = 0.0;

    for (size_t i = 0; i < 5000; i++)
    {

    t1 = std::chrono::high_resolution_clock::now();
    result = new CySolverResult(num_y, num_y, expected_size);

    solver = new RK45(
        diffeq_func,
        t_start,
        t_end,
        y0_ptr,
        num_y,
        false,
        0,
        NULL,
        0,
        2000,
        rtol,
        atol,
        NULL,
        NULL,
        INF,
        0.0
    );
    step_success = false;
    success = false;
    
    solver->reset();
    result->save_data(solver->t_now, solver->y_now_ptr, solver->dy_now_ptr, true, solver->error_code, solver->message_ptr);
        
    while ((solver->status == 0) and (solver->error_code == 0))
    { // and (result.cyres_ptr.error_code == 0) 
        solver->take_step();
        if (solver->status < 0)
        {
            break;
        }

        step_success = (solver->status == 0) or (solver->status == 1);
        if (step_success)
        {
            result->save_data(solver->t_now, solver->y_now_ptr, solver->dy_now_ptr, step_success, solver->error_code, solver->message_ptr);
        }
        else {
            break;
        }
        
        if (solver->status == 1)
        {
            break;
        }
    }

    if (solver->status == 1)
    {
        success = true;
    }

    final_size = result->size;
    sol_size = final_size * 2;

    delete result;
    delete solver;

    t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::micro> ms_double = t2 - t1;
    running_sum += ms_double.count();
    //std::cout << ms_double.count() << "ms\n";
    }

    std::cout << "\nSIZE: " << final_size << std::endl;
    std::cout << "\n\n AVERAGE: " << running_sum / 5000 << "ms\n" << std::endl;
    //break;
    /*std::cout << "Done! Final size: " << final_size << std::endl;

    std::ofstream datastream;

    datastream.open("out.dat");
    if (!datastream)
    {
        std::cout << "Could not open data file." << std::endl;
    }
    datastream.precision(6);

    for (size_t i = 0; i < final_size; i++)
    {
        datastream << result->time_domain_ptr[0][i] << ", " << result->solution_ptr[0][2 * i] << ", " << result->solution_ptr[0][2 * i + 1] << std::endl;
    }
    datastream.close();
    std::cout << "Data saved." << std::endl;*/

    //std::cin.get();


    }
    return 0;
}