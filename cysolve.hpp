#pragma once

#include <memory>

#include "common.hpp"
#include "rk.hpp"
#include "cysolver.hpp"

template <typename T>
void find_cysolver_and_solve(
    DiffeqFuncType diffeq_ptr,
    std::shared_ptr<CySolverResult> solution_ptr,
    const double t_start,
    const double t_end,
    const double* y0_ptr,
    const unsigned int num_y,
    // General optional arguments
    const unsigned int num_extra,
    const double* args_ptr,
    // rk optional arguments
    const size_t max_num_steps,
    const size_t max_ram_MB,
    const double rtol,
    const double atol,
    const double* rtols_ptr,
    const double* atols_ptr,
    const double max_step_size,
    const double first_step_size
);

template <typename T>
void find_pysolver_and_solve(
    // Cython class instance used for pyhook
    PyObject* cython_extension_class_instance,
    std::shared_ptr<CySolverResult> solution_ptr,
    const double t_start,
    const double t_end,
    const double* y0_ptr,
    const unsigned int num_y,
    // General optional arguments
    const unsigned int num_extra,
    const double* args_ptr,
    // rk optional arguments
    const size_t max_num_steps,
    const size_t max_ram_MB,
    const double rtol,
    const double atol,
    const double* rtols_ptr,
    const double* atols_ptr,
    const double max_step_size,
    const double first_step_size
);

std::shared_ptr<CySolverResult> cysolve_ivp(
    DiffeqFuncType diffeq_ptr,
    const double* t_span_ptr,
    const double* y0_ptr,
    const unsigned int num_y,
    const unsigned int method,
    // General optional arguments
    const size_t expected_size = 0,
    const unsigned int num_extra = 0,
    const double* args_ptr = nullptr,
    // rk optional arguments
    const size_t max_num_steps = 0,
    const size_t max_ram_MB = 2000,
    const double rtol = 1.0e-3,
    const double atol = 1.0e-6,
    const double* rtols_ptr = nullptr,
    const double* atols_ptr = nullptr,
    const double max_step_size = MAX_STEP,
    const double first_step_size = 0.0
);