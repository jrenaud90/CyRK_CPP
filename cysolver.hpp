#pragma once

#include <cstring>

#include <memory>

#include <Python.h>

#include "common.hpp"
#include "cysolution.hpp"

class CySolverBase {


// Attributes
protected:
    // ** Attributes **

    // Time variables
    double t_start       = 0.0;
    double t_end         = 0.0;
    double t_old         = 0.0;
    double t_delta       = 0.0;
    double t_delta_abs   = 0.0;
    double direction_inf = 0.0;

    // Dependent variables
    unsigned int num_dy = 0;
    double num_y_dbl  = 0.0;
    double num_y_sqrt = 0.0;
    // The size of the stack allocated tracking arrays is equal to the maximum allowed `num_y` (25).
    double y0[Y_LIMIT]    = { 0.0 };
    double y_old[Y_LIMIT] = { 0.0 };
    double y_now[Y_LIMIT] = { 0.0 };
    // For dy, both the dy/dt and any extra outputs are stored. So the maximum size is `num_y` (25) + `num_extra` (25)
    double dy_old[DY_LIMIT] = { 0.0 };
    double dy_now[DY_LIMIT] = { 0.0 };

    // dy_now_ptr and y_now_ptr are declared in public.
    double* y0_ptr     = &y0[0];
    double* y_old_ptr  = &y_old[0];
    double* dy_old_ptr = &dy_old[0];

    // Integration step information
    size_t max_num_steps = 0;

    // Differential equation information
    const double* args_ptr    = nullptr;
    DiffeqFuncType diffeq_ptr = nullptr;

    // Information on capturing extra information during integration.
    int num_extra = 0;

    // Keep bools together to reduce size
    bool direction_flag = false;
    bool reset_called   = false;
    bool capture_extra  = false;
    bool user_provided_max_num_steps = false;

public:
    // Status attributes
    int status = -999;

    // Meta data
    unsigned int num_y = 0;

    // Result storage
    std::shared_ptr<CySolverResult> storage_ptr = nullptr;

    // State attributes
    size_t len_t = 0;
    double t_now = 0.0;
    double* y_now_ptr = &y_now[0];
    double* dy_now_ptr = &dy_now[0];

    // PySolver Attributes
    bool use_pysolver = false;
    PyObject* cython_extension_class_instance = nullptr;


// Methods
protected:
    virtual void p_step_implementation();

public:
    CySolverBase();
    virtual ~CySolverBase();
    CySolverBase(
        // Input variables
        DiffeqFuncType diffeq_ptr,
        std::shared_ptr<CySolverResult> const storage_ptr,
        const double t_start,
        const double t_end,
        const double* const y0_ptr,
        const unsigned int num_y,
        const unsigned int num_extra = 0,
        const double* const args_ptr = nullptr,
        const size_t max_num_steps = 0,
        const size_t max_ram_MB = 2000
    );
    
    bool check_status() const;
    virtual void reset();
    void diffeq();
    void take_step();
    void change_storage(std::shared_ptr<CySolverResult> new_storage_ptr, bool auto_reset = true);

    // PySolver methods
    void set_cython_extension_instance(PyObject* cython_extension_class_instance);
    void py_diffeq();
};
