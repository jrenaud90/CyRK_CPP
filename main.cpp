/* This file and its functions are just used for crude testing. It is not considered truely apart of CyRK's C++ backend. */
#include <iostream>
#include <fstream>
#include <chrono>
#include <cstring>
#include <vector>
# define M_PI           3.14159265358979323846 

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

static void large_numy_diffeq(double* dy_ptr, double t, double* y_ptr, char* args_ptr, PreEvalFunc pre_eval_func)
{
    size_t num_y = 10000;

    // Fix: reinterpret_cast is required for pointer-to-object conversions between unrelated types
    double* args_dbl_ptr = reinterpret_cast<double*>(args_ptr);
    double decay_rate = args_dbl_ptr[0];
    double forcing_scale = args_dbl_ptr[1];

    // This diffeq converges so should be stable
    for (size_t i = 0; i < num_y; i++)
    {
        decay_rate *= 0.9999;

        if (i < (num_y - 1))
        {
            dy_ptr[i] = decay_rate * y_ptr[i] * sin(2 * M_PI * t / 5.0 + y_ptr[i + 1] / 50.0);
        }
        else
        {
            dy_ptr[i] = decay_rate * y_ptr[i] * sin(2 * M_PI * t / 5.0);
        }
    }
}

static void large_numy_simple_diffeq(double* dy_ptr, double t, double* y_ptr, char* args_ptr, PreEvalFunc pre_eval_func)
{
    const size_t num_y = 10000;

    dy_ptr[0] = sin(2 * M_PI * t / 10.0);

    // Set all else to zero.
    memset(&dy_ptr[0], 0, num_y - 1);
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
    //DiffeqFuncType diffeq_func = large_numy_diffeq;
    DiffeqFuncType diffeq_func = large_numy_simple_diffeq;
    if (num_extra > 0)
    {
        DiffeqFuncType diffeq_func = test_extra_diffeq;
    }

    double t_start = 0.0;
    int max_i = 300;

    //std::vector<double> y0_vec = std::vector<double>(2);
    //y0_vec[0] = 20.0; // Initial condition for y0
    //y0_vec[1] = 20.0; // Initial condition for y0
    
    std::vector<double> y0_vec = std::vector<double>(10000);
    for (size_t i = 0; i < 10000; i++)
    {
        y0_vec[i] = 100.0;
    }

    std::string msg = "No message";

    unsigned int num_y = 2;

    std::vector<double> rtols_vec = std::vector<double>(1);
    rtols_vec[0] = 1.0e-7; // Relative tolerance for y0 and y1
    std::vector<double> atols_vec = std::vector<double>(1);
    atols_vec[0] = 1.0e-8; // Absolute tolerance for y0 and y1

    size_t expected_size = 0;

    double running_sum = 0.0;
    double running_sum_ps = 0.0;
    size_t final_size;
    size_t sol_size;
    size_t num_interps = 0;
    std::vector<char> args_vec = std::vector<char>(2 * sizeof(double));
    double* args_dbl_ptr = reinterpret_cast<double*>(args_vec.data());
    args_dbl_ptr[0] = -0.5;
    args_dbl_ptr[1] = 1.0e-5;

    std::vector<Event> events_vec = std::vector<Event>(0);

    std::vector<double> t_eval_vec = std::vector<double>(0);
    if (len_t_eval > 0)
    {
        t_eval_vec = linspace(t_start, t_end, len_t_eval);
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    auto t2 = std::chrono::high_resolution_clock::now();
    int k = 0;
    int k_max = 20; // 15
    double total_runner = 0.0;
    double total_runner_ps = 0.0;
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
        running_sum_ps = 0.0;

        double cheat_check = t_end; // uncomment if you want to make sure things are changing every step t_end + 100 * k;

        for (size_t i = 0; i < max_i; i++)
        {
            if (use_resets)
            {
                result_ptr = result_uptr_forresets.get();
                t1 = std::chrono::high_resolution_clock::now();
                baseline_cysolve_ivp_noreturn(
                    result_ptr,
                    diffeq_func,
                    t_start,
                    cheat_check,
                    y0_vec,
                    expected_size,
                    num_extra,
                    args_vec, // args_vec
                    100000,
                    2000,
                    save_dense,
                    t_eval_vec,
                    nullptr, // pre_eval_func
                    events_vec,
                    rtols_vec,
                    atols_vec,
                    INF, // max_step_size
                    0.0, // first_step_size
                    true  // Force retain solver.
                );
                t2 = std::chrono::high_resolution_clock::now();
            }
            else
            {
                t1 = std::chrono::high_resolution_clock::now();
                result_uptr_regular = baseline_cysolve_ivp(
                    diffeq_func,
                    t_start,
                    t_end,
                    y0_vec,
                    method,
                    expected_size,
                    num_extra,
                    args_vec,
                    100000,
                    2000,
                    save_dense,
                    t_eval_vec,
                    nullptr, // pre_eval_func
                    events_vec,
                    rtols_vec,
                    atols_vec,
                    INF, // max_step_size
                    0.0, // first_step_size
                    true  // Force retain solver.
                );
                t2 = std::chrono::high_resolution_clock::now();
                result_ptr = result_uptr_regular.get();
            }

            final_size = result_ptr->size;
            msg = result_ptr->message;
            /*if (result_ptr->success)
            {
                num_interps = result_ptr->num_interpolates;
                t_at_20 = result_ptr->time_domain_vec[4];
                y_at_20_ptr = &result_ptr->solution[(num_y + num_extra) * 4];
                result_ptr->call(t_at_20, y_int_at_20_ptr);
            }
            else
            {
                num_interps = -1;
                t_at_20 = NULL;
                y_at_20_ptr = nullptr;
                y_int_at_20_ptr[0] = NULL;
                y_int_at_20_ptr[1] = NULL;
            }*/

            //std::cout << result->message_ptr << std::endl;
            sol_size = final_size * 2;

            std::chrono::duration<double, std::micro> ms_double = t2 - t1;
            running_sum += ms_double.count();
            running_sum_ps += ms_double.count() / final_size;
            //std::cout << ms_double.count() << "ms\n";
        }

        std::cout << "SIZE: " << final_size << std::endl;
        std::cout << "NUM INTERPS: " << num_interps << std::endl;
        /*if (y_at_20_ptr)
        {
            std::cout << "20 t = " << t_at_20 << "; y0 = " << y_at_20_ptr[0] << ", y1 = " << y_at_20_ptr[1] << std::endl;
        }
        else
        {
            std::cout << "20 t = " << t_at_20 << "; y0 = " << NULL << ", y1 = " << NULL << std::endl;
        }
        std::cout << "DENSE INTERP:: " << "y0 = " << y_int_at_20_ptr[0] << ", y1 = " << y_int_at_20_ptr[1] << std::endl;*/
        std::cout << "Message: " << msg << std::endl;
        std::cout << "AVERAGE: " << running_sum / max_i << "us\n" << std::endl;
        std::cout << "AVERAGE PS: " << running_sum_ps / max_i << "us\n" << std::endl;
        if (k > 3)
        {
            total_runner += running_sum / max_i;
            total_runner_ps += running_sum_ps / max_i;
        }
        k += 1;
        //break;
    }

    std::cout << "Done! Final Avg: " << total_runner / (k_max - 4) << std::endl;
    std::cout << "      Final Avg PS: " << total_runner_ps / (k_max - 4) << std::endl << "Saving data..." << std::endl;


    std::ofstream datastream;

    datastream.open("out.dat");
    if (!datastream)
    {
        std::cout << "Could not open data file." << std::endl;
    }
    datastream.precision(32);

    /*if (result_ptr)
    {
        for (size_t i = 0; i < final_size; i++)
        {
            datastream << result_ptr->time_domain_vec[i] << ", " << result_ptr->solution[2 * i] << ", " << result_ptr->solution[2 * i + 1] << std::endl;
        }
    }*/
    datastream.close();
    std::cout << "Data saved." << std::endl;

    //std::cin.get();
}

int main(){
    /* This file and its functions are just used for crude testing. It is not considered truely apart of CyRK's C++ backend. */

    /* Speed Runs (JPR Work Desktop):
    * All times in microseconds.
    * No Dense; No t_eval; t_end = 500.0; method=RK45
    *   Steps: 3535; t20 = 0.6101; y0 = 33.18; y1 = 14.91
    *   Older:   755.8; 610.9
    *   v0.13.0: 932.7; 933.7; 934.3
    *   v0.14.0 (no resets): 373; 372.2
    *   v0.14.0 (resets):    375.4; 372.95
    *   v0.15.0 (resets):    372.6; 371.6
    * 
    * No Dense; No t_eval; t_end = 500.0; method=DOP853
    *   Steps: 1236; t20 = 1.866; y0 = 96.679; y1 = 18.98
    *   v0.14.0 (resets):    284.7; 284.914
    *   v0.15.0 (resets):    281.9; 282.7
    *
    * Short time step; t_end = 5.0; method=RK45
    *   v0.14.0 (no resets): 6.1; 6.1; 5.7
    *   v0.14.0 (resets):    4.87; 5.01; 4.95
    *   v0.15.0 (resets):    4.75; 4.59; 4.66
    * 
    * With Dense; No t_eval
    *   Older:   1089.3; 1094.1
    *   v0.13.0: 1812.1; 1960.9; 1847.1
    *   v0.14.0 (no resets): 1937.7; 1906.1
    *   v0.14.0 (resets):    1546.8; 1578.4
    *
    * With Dense; With 2x t_eval
    *   Older:   1246.6; 1239.9
    *   v0.13.0: 2496.4; 2468.3
    *   v0.14.0 (no resets): 2626.7; 2622.0
    *   v0.14.0 (resets):    2449.7; 2476.2
    * 
    *
    * With Dense; With 0.5x t_eval
    *   Older:   1135.5; 1136.5
    *   v0.13.0: 2170.5; 2198.8
    *   v0.14.0 (no resets): 2267.1; 2256.0
    *   v0.14.0 (resets):    1965.9; 1995.0
    *
    * No Dense; With 0.5x t_eval
    *   Older:   822.6; 820.4
    *   v0.13.0: 1336.3; 1303.9
    *   v0.14.0 (no resets): 1224.0; 1220.5
    *   v0.14.0 (resets):    1249.6; 1249.2
    *
    * No Dense; With 2x t_eval
    *   Older:   903; 901.4
    *   v0.13.0: 1631; 1631
    *   v0.14.0 (no resets): 1595.0; 1599.4
    *   v0.14.0 (resets):    1577.0; 1569.4
    *
    * No Dense; extra on
    *   v0.14.0 (no resets): 934.5; 939.0
    *   v0.14.0 (resets):    940.5; 938.8
    * 
    * Dense On; extra on
    *   v0.14.0 (no resets): 1961.8; 2068.4
    *   v0.14.0 (resets):    1872.0; 1611.8
    * 
    * No Dense; With 2x t_eval; With extra on
    *   Older:   1000.79us; 952.867; 941.13
    *   v0.13.0: broken?
    *   v0.14.0 (no resets): 1866.6; 1877.1
    *   v0.14.0 (resets):    1825.7; 1820.6
    *
    */

    test_regular(
        50.0, // t_end
        false, // Dense
        0,     // len t_eval 7070 == 2x; 1767 == 0.5x for tspan of (0, 500))
        0,     // num extra
        ODEMethod::DOP853,      // Method
        true  // use_resets
    );

    // 1333.07

    // Large num_y (simple) tests 2026-02-18
    // t_end = 50.0; dense = false; num_extra = 0; len_t_eval = 0; RK45; resets = True
    // Size = 86
    // v1.0 - Avg 13782.6   PS: 160.262
    // v1.1 - Avg 13120.9   PS: 152.569
    // v1.2 - Avg 18667.4   PS: 217.063
    // v1.3 - Avg 12075.1   PS: 140.408
    // vC2.0 (new diffeq) - Avg 71664.3   PS: 421.555
    // vS2.0 (new diffeq) - Avg 925.595   PS: 92.5595  ; 935.196 PS: 93.5196
    // vS2.1 - Avg 925.595   PS: 92.5595
    // vS2.1 - Avg 905.809   PS: 90.5809
    // vS2.2 - Avg 892.034   PS: 89.2034
    // vS2.3 - Avg 934.32   PS: 93.432
    // vS2.3 - Avg 881.91   PS: 88.191
    // vS-DEN2.3 - Avg 1679.82   PS: 167.982
    // vS-DEN2.3 - Avg 1626.39   PS: 162.639
    // vS-DEN2.4 - Avg 1507.08   PS:  150.708
    // vS-DOP2.4 - Avg 2439.71   PS:  243.971
    // vS-DOP-DEN2.5 - Avg 4909.61   PS:  490.961
    // vS-DOP2.4 - Avg 2497.13   PS:  249.713

    // vS-RK4.0 - Avg  909.039   PS:  90.9039
    // vS-DOP4.0 - Avg 2411.5   PS:   241.15
    // vS-RK4.1 - Avg  723.064   PS:  72.3064
    // vS-DOP4.1 - Avg 1764.29   PS:   176.429
    // vS-RK4.2 - Avg  651.648   PS:   65.1648
    // vS-DOP4.2 - Avg 1643.67   PS:   164.367
    // vS-RK-DENSE4.2 - Avg  1290.7   PS:   129.07
    // vS-DOP-DENSE4.2 - Avg  4069.5   PS:   406.95
    // vS-RK-DENSE4.3 - Avg  1090.9   PS:   109.09
    // vS-DOP-DENSE4.3 - Avg   3449.36   PS:   344.936

    // vS-RK-DENSE4.4 - Avg  1141.5   PS:  114.15
    // vS-DOP-DENSE4.4 - Avg   3434.51   PS:   343.451
    // vS-RK4.4 - Avg  674.137   PS:   67.4137
    // vS-DOP4.4 - Avg  1667.48   PS:   166.748
    // vS-RK4.5 - Avg  543.365   PS:   54.3365
    // vS-DOP4.5 - Avg  1547.58   PS:   154.758
    
    // CyRK 0.17.0
    // RK17.0a        - Avg  546.148    PS:   54.6148
    // DOP17.0a       - Avg  1562       PS:   156.2
    // RK-DENSE17.0a  - Avg  1019.03    PS:  101.903
    // DOP-DENSE17.0a - Avg  3369.33    PS:   336.933

    return 0;
}