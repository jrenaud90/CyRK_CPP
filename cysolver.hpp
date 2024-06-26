#pragma once

#include <vector>
#include <limits>
#include <cmath>
#include <cstring>
#include <cstdio>

#include "rk_constants.h"

// Integration Constants
// Multiply steps computed from asymptotic behaviour of errors by this.
const double SAFETY = 0.9;  // Error coefficient factor (1 == use calculated error; < 1 means be conservative).
const double MIN_FACTOR = 0.2;  // Minimum allowed decrease in a step size.
const double MAX_FACTOR = 10.;  // Maximum allowed increase in a step size.
const double INF = std::numeric_limits<double>::infinity();
const double MAX_STEP = INF;
const double EPS = std::numeric_limits<double>::epsilon();
const double EPS_10 = EPS * 10.0;
const double EPS_100 = EPS * 100.0;
const size_t MAX_SIZET_SIZE = std::numeric_limits<size_t>::max();
const size_t MAX_INT_SIZE = std::numeric_limits<int>::max();

// Solution parameters
const double DYNAMIC_GROWTH_RATE = 1.618;
// To learn why the golden ratio is used, read this:
// https://stackoverflow.com/questions/1100311/what-is-the-ideal-growth-rate-for-a-dynamically-allocated-array
const double SIZE_MAX_DBL = 0.99 * SIZE_MAX;


void find_max_num_steps(
    size_t num_y,
    size_t num_extra,
    size_t max_num_steps,
    size_t max_ram_MB,
    bool capture_extra,
    bool* user_provided_max_num_steps,
    size_t* max_num_steps_touse);


typedef void (*DiffeqFuncType)(double*, double, double*, double*);

class CySolverResult {

protected:
    // ** Attributes **
    // Message storage
    char message[256];

    // Metadata
    size_t num_y;
    size_t num_dy;
    double num_dy_dbl;

    // Current storage information
    size_t storage_capacity;

    // ** Methods **
    void expand_storage();
public:
    // ** Attributes **
    // Status information
    char* message_ptr;
    bool success;
    int error_code;
    size_t size;

    // Storage for arrays
    std::vector<double> time_domain;
    std::vector<double> solution;
    std::vector<double>* time_domain_ptr;
    std::vector<double>* solution_ptr;

    // ** Methods **
    CySolverResult();
    CySolverResult(size_t num_y, size_t num_dy, size_t expected_size);
    void save_data(double new_t, double* new_solution_y, double* new_solution_dy, bool success, int error_code, char* message);
};

class CySolverBase {

protected:
    // ** Attributes **
    // Status variables
    bool reset_called;
    char message[512];

    // Time variables
    double t_start;
    double t_end;
    double t_old;
    double t_delta;
    double t_delta_abs;
    double direction_inf;
    bool direction_flag;

    // Dependent variables
    size_t num_y;
    double num_y_dbl;
    double num_y_sqrt;
    size_t num_dy;
    // The size of the stack allocated tracking arrays is equal to the maximum allowed `num_y` (50).
    double y0[50] = { std::nan("") };
    double y_old[50] = { std::nan("") };
    double y_now[50] = { std::nan("") };
    // For dy, both the dy/dt and any extra outputs are stored. So the maximum size is `num_y` (50) + `num_extra` (50)
    double dy_old[100] = { std::nan("") };
    double dy_now[100] = { std::nan("") };

    // dy_now_ptr and y_now_ptr are declared in public.
    double* y0_ptr = &y0[0];
    double* y_old_ptr = &y_old[0];
    double* dy_old_ptr = &dy_old[0];

    // Integration step information
    size_t len_t;
    size_t max_num_steps;
    bool user_provided_max_num_steps;

    // Information on capturing extra information during integration.
    bool capture_extra;
    size_t num_extra;

    // Differential equation information
    double* args_ptr;
    DiffeqFuncType diffeq_ptr;

    // ** Methods **
    virtual void step_implementation();

public:

    // ** Attributes **
    int status;
    // Error Codes:
    // -1 : Attribute Error
    // -2 : Value Error
    // -3 : Not Implemented Error
    // -4 : <Reserved>
    // -5 : Runtime error
    int error_code;
    char* message_ptr;

    // State properties
    double t_now;
    double* y_now_ptr = &y_now[0];
    double* dy_now_ptr = &dy_now[0];

    // ** Methods **
    // Constructors
    CySolverBase();
    CySolverBase(
        // Input variables
        DiffeqFuncType diffeq_ptr,
        double t_start,
        double t_end,
        double* y0_ptr,
        size_t num_y,
        bool capture_extra = false,
        size_t num_extra = 0,
        double* args_ptr = nullptr,
        size_t max_num_steps = 0,
        size_t max_ram_MB = 2000
    );
    void diffeq();
    virtual void reset();
    void take_step();
};


class RKSolver : public CySolverBase {

protected:
    // ** Attributes **

    // Step globals
    const double error_safety = SAFETY;
    const double min_step_factor = MIN_FACTOR;
    const double max_step_factor = MAX_FACTOR;

    // RK constants
    int order;
    int error_estimator_order;
    double error_exponent;
    size_t n_stages;
    size_t len_Acols;
    size_t len_C;
    size_t nstages_numy;

    // Pointers to RK constant arrays
    const double* C_ptr;
    const double* A_ptr;
    const double* B_ptr;
    const double* E_ptr;
    const double* E3_ptr;
    const double* E5_ptr;
    const double* P_ptr;
    const double* D_ptr;
    double* K_ptr;

    // K is not const. Its values are stored in an array that is held by this class.
    double K[1] = { std::nan("") };

    // Tolerances
    // For the same reason num_y is limited, the total number of tolerances are limited.
    double rtols[50] = { std::nan("") };
    double atols[50] = { std::nan("") };
    double* rtols_ptr = &rtols[0];
    double* atols_ptr = &atols[0];
    bool use_array_rtols;
    bool use_array_atols;

    // Step size parameters
    double step;
    double step_size_old;
    double step_size;
    double max_step_size;
    double user_provided_first_step_size;

    // Error estimate
    double error_norm;
    void estimate_error();
    void step_implementation() override;

public:
    RKSolver();
    RKSolver(
        // Input variables
        DiffeqFuncType diffeq_ptr,
        double t_start,
        double t_end,
        double* y0_ptr,
        size_t num_y,
        bool capture_extra = false,
        size_t num_extra = 0,
        double* args_ptr = nullptr,
        size_t max_num_steps = 0,
        size_t max_ram_MB = 2000,
        double rtol = 1.0e-3,
        double atol = 1.0e-6,
        double* rtols_ptr = nullptr,
        double* atols_ptr = nullptr,
        double max_step_size = MAX_STEP,
        double first_step_size = 0.0
    );
    virtual void reset() override;
    void calc_first_step_size();
};


class RK23 : public RKSolver {

protected:
    double K[4 * 50] = { 0.0 };

public:
    // Copy over base class constructors
    using RKSolver::RKSolver;

    virtual void reset() override;
};


class RK45 : public RKSolver {

protected:
    double K[7 * 50] = { 0.0 };

public:
    // Copy over base class constructors
    using RKSolver::RKSolver;
    virtual void reset() override;
};