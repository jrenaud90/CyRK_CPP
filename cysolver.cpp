#include "cysolver.hpp"

// !!!
// Uncomment these dummy methods if working outside of CyRK and you just want the program to compile and run for testing/developing the C++ only code.

bool import_CyRK__cy__pysolver_cyhook()
{
    return true;
}

int call_diffeq_from_cython(PyObject* x, DiffeqMethod y)
{
    return 1;
}

void Py_XINCREF(PyObject* x)
{
}

void Py_XDECREF(PyObject* x)
{
}

// Constructors
CySolverBase::CySolverBase() {}
CySolverBase::CySolverBase(
    DiffeqFuncType diffeq_ptr,
    std::shared_ptr<CySolverResult> storage_ptr,
    const double t_start,
    const double t_end,
    const double* const y0_ptr,
    const unsigned int num_y,
    const unsigned int num_extra,
    const double* const args_ptr,
    const size_t max_num_steps,
    const size_t max_ram_MB,
    const bool use_dense_output,
    const double* t_eval,
    const size_t len_t_eval) :
    status(0),
    num_y(num_y),
    num_extra(num_extra),
    t_start(t_start),
    t_end(t_end),
    storage_ptr(storage_ptr),
    diffeq_ptr(diffeq_ptr),
    args_ptr(args_ptr),
    use_dense_output(use_dense_output),
    t_eval_ptr(t_eval),
    len_t_eval(len_t_eval)
{
    // Parse inputs
    this->capture_extra = num_extra > 0;

    // Setup storage
    this->storage_ptr->update_message("CySolverBase Initializing.");

    // Check for errors
    if (this->num_extra > (DY_LIMIT - Y_LIMIT))
    {
        this->status = -9;
        this->storage_ptr->error_code = -1;
        this->storage_ptr->update_message("CySolverBase Attribute Error: `num_extra` exceeds the maximum supported size.");
    }

    if (this->num_y > Y_LIMIT)
    {
        this->status = -9;
        this->storage_ptr->error_code = -1;
        this->storage_ptr->update_message("CySolverBase Attribute Error: `num_y` exceeds the maximum supported size.");
    }
    else if (this->num_y == 0)
    {
        this->status = -9;
        this->storage_ptr->error_code = -1;
        this->storage_ptr->update_message("CySolverBase Attribute Error: `num_y` = 0 so nothing to integrate.");
    }

    // Parse y values
    this->num_y_dbl  = (double)this->num_y;
    this->num_y_sqrt = std::sqrt(this->num_y_dbl);
    this->num_dy     = this->num_y + this->num_extra;
    // Make a copy of y0
    std::memcpy(this->y0_ptr, y0_ptr, sizeof(double) * this->num_y);

    // Parse time information
    this->t_delta = t_end - t_start;
    this->t_delta_abs = std::fabs(this->t_delta);
    if (this->t_delta >= 0.0)
    {
        // Forward integration
        this->direction_flag = true;
        this->direction_inf = INF;
    }
    else {
        // Backward integration
        this->direction_flag = false;
        this->direction_inf = -INF;
    }

    // Parse maximum number of steps
    MaxNumStepsOutput max_num_steps_output = find_max_num_steps(
        this->num_y,
        num_extra,
        max_num_steps,
        max_ram_MB
    );
    this->user_provided_max_num_steps = max_num_steps_output.user_provided_max_num_steps;
    this->max_num_steps               = max_num_steps_output.max_num_steps;

    // Bind diffeq to C++ version
    this->diffeq = &CySolverBase::cy_diffeq;

    // Parse t_eval
    if (this->t_eval_ptr && this->len_t_eval > 0)
    {
        this->use_t_eval = true;
        if (!this->direction_flag)
        {
            // TODO: add in support for backwards integration when t_eval is provided.
            this->status = -70;
            this->storage_ptr->error_code = -70;
            this->storage_ptr->update_message("Not Implemented Error: Can only do forward integration when t_eval is provided.");
        }
    }
}


// Destructors
CySolverBase::~CySolverBase()
{
    this->storage_ptr = nullptr;
    if (this->use_pysolver)
    {
        // Decrease reference count on the cython extension class instance
        Py_XDECREF(this->cython_extension_class_instance);
    }
}


// Protected methods
void CySolverBase::p_estimate_error()
{
    // Overwritten by subclasses.
}

void CySolverBase::p_step_implementation()
{
    // Overwritten by subclasses.
}

// Public methods
void CySolverBase::calc_first_step_size()
{
    // Overwritten by subclasses.
}

bool CySolverBase::check_status() const
{
    // If the solver is not in state 0 then that is an automatic rejection.
    if (this->status != 0)
    {
        return false;
    }

    // Otherwise, check if the solution storage is in an error state.
    if (this->storage_ptr) [[likely]]
    {
        if (this->storage_ptr->error_code != 0)
        {
            return false;
        }
    }

    // If we reach here then we should be good to go.
    return true;
}

void CySolverBase::cy_diffeq() noexcept
{
    // Call c function
    this->diffeq_ptr(this->dy_now_ptr, this->t_now_ptr[0], this->y_now_ptr, this->args_ptr);
}

void CySolverBase::reset()
{
    this->status = 0;
    this->reset_called = false;

    // Reset time
    this->t_now_ptr[0] = this->t_start;
    this->t_old = this->t_start;
    this->len_t = 1;

    // Reset ys
    std::memcpy(this->y_now_ptr, this->y0_ptr, sizeof(double) * this->num_y);
    std::memcpy(this->y_old_ptr, this->y0_ptr, sizeof(double) * this->num_y);

    // Call differential equation to set dy0
    this->diffeq(this);

    // Update dys
    std::memcpy(this->dy_old_ptr, this->dy_now_ptr, sizeof(double) * this->num_y);

    // Initialize storage
    this->storage_ptr->reset();
    this->storage_ptr->update_message("CySolverStorage reset, ready for data.");

    // If t_eval is set then don't save initial conditions. They will be captured during stepping.
    if (!this->use_t_eval)
    {
        // Store initial conditions
        this->storage_ptr->save_data(this->t_now_ptr[0], this->y_now_ptr, this->dy_now_ptr);
    }
    
    // Construct interpolator using t0 and y0 as its data point
    if (this->use_dense_output)
    {
        CySolverDense* dense_output_ptr = this->p_dense_output_heap();
        // Save interpolator
        this->storage_ptr->save_dense(this->t_now_ptr[0], dense_output_ptr);
    }

    // Prep for t_eval
    this->t_eval_index_old = 0;

    // Done with reset
    this->reset_called = true;
}

#include <cstdio>

void CySolverBase::take_step()
{    
    if (!this->reset_called) [[unlikely]]
    {
        // Reset must be called first.
        this->reset();
    }

    bool skip_t_eval = false;

    if (!this->status)
    {
        if (this->t_now_ptr[0] == this->t_end) [[unlikely]]
        {
            // Integration finished
            this->t_old = this->t_end;
            this->status = 1;
        }
        else if (this->len_t >= this->max_num_steps) [[unlikely]]
        {
            if (this->user_provided_max_num_steps)
            {
                // Maximum number of steps reached (as set by user).
                this->status = -2;
            }
            else {
                // Maximum number of steps reached (as set by RAM limitations).
                this->status = -3;
            }
        }
        else [[likely]]
        {
            // ** Make call to solver's step implementation **
            bool save_data = true;
            bool prepare_for_next_step = true;

            this->p_step_implementation();
            this->len_t++;

            // Take care of dense output and t_eval
            if (this->use_dense_output)
            {
                // We need to save many dense interpolators to storage. So let's heap allocate them.
                CySolverDense* dense_output_heap_ptr = this->p_dense_output_heap();
                // Save to storage array.
                this->storage_ptr->save_dense(this->t_now_ptr[0], dense_output_heap_ptr);
            }

            if (this->use_t_eval && !skip_t_eval)
            {
                // Don't save data at the end
                save_data = false;

                // We are not saving interpolators to storage but we still need one to work on t_eval. 
                // We will only ever need 1 interpolator per step. So let's just stack allocate that one.
                CySolverDense dense_output(
                    this->integration_method,
                    this->t_old,
                    this->t_now_ptr[0],
                    this->y_old_ptr,
                    this->num_y,
                    0 // Fake Q order just for consistent constructor call
                    );
                // Update the dense output class with integrator-specific data
                this->p_dense_output_stack(dense_output);

                // Need to step through t_eval and call dense to determine correct data at each t_eval step.
                // Find the first index in t_eval that is close to current time.
                size_t t_eval_index_new = 1 + binary_search_with_guess(this->t_now_ptr[0], this->t_eval_ptr, this->len_t_eval, this->t_eval_index_old);
                if (t_eval_index_new >= this->len_t_eval)
                {
                    t_eval_index_new = this->len_t_eval;
                    // We are done with t_eval. Skip it from now on. 
                    skip_t_eval = true;
                }                
                
                // Check if there are any t_eval steps between this new index and the last index.
                int t_eval_index_delta = (int)t_eval_index_new - (int)this->t_eval_index_old;
                // If t_eval_index_delta == 0 then there are no new interpolations required between the last integration step and now.
                // ^ In this case do not save any data, we are done with this step.

                if (t_eval_index_delta > 0)
                {
                    // There are steps we need to interpolate over.
                    // Start with the old time and add t_eval step sizes until we are done.
                    // Create a y array and dy_array to use during interpolation
                    double y_interp[Y_LIMIT] = { };
                    double* y_interp_ptr     = &y_interp[0];

                    // If capture extra is set to true then we need to hold onto a copy of the current state
                    // The current state pointers must be overwritten if extra output is to be captured.
                    // However we need a copy of the current state pointers at the end of this step anyways. So just
                    // store them now and skip storing them later.

                    if (this->capture_extra)
                    {
                        // We need to copy the current state of y, dy, and time
                        this->t_old = this->t_now_ptr[0];
                        std::memcpy(this->y_old_ptr, this->y_now_ptr, sizeof(double) * this->num_y);
                        std::memcpy(this->dy_old_ptr, this->dy_now_ptr, sizeof(double) * this->num_dy);

                        // Don't update these again at the end
                        prepare_for_next_step = false;
                    }

                    for (size_t i = 0; i < t_eval_index_delta; i++)
                    {
                        double t_interp = this->t_eval_ptr[this->t_eval_index_old + i];

                        // Call the interpolator using this new time value.
                        dense_output.call(t_interp, y_interp_ptr);

                        if (this->capture_extra)
                        {
                            // If the user want to capture extra output then we also have to call the differential equation to get that extra output.
                            // To do this we need to hack the current integrators t_now, y_now, and dy_now.
                            // TODO: This could be more efficient if we just changed pointers but since the PySolver only stores y_now_ptr, dy_now_ptr, etc at initialization, it won't be able to see changes to new pointer. 
                            // So for now we have to do a lot of copying of data.

                            // Copy the interpreted y onto the current y_now_ptr. Also update t_now
                            this->t_now_ptr[0] = t_interp;
                            std::memcpy(this->y_now_ptr, y_interp_ptr, sizeof(double) * this->num_y);

                            // Call diffeq to update dy_now_ptr with the extra output.
                            this->diffeq(this);
                        }
                        // Save interpolated data to storage. If capture extra is true then dy_now holds those extra values. If it is false then it won't hurt to pass dy_now to storage.
                        this->storage_ptr->save_data(t_interp, y_interp_ptr, this->dy_now_ptr);
                    }
                }
                // Update the old index for the next step
                this->t_eval_index_old = t_eval_index_new;
            }
            if (save_data)
            {
                // No data has been saved from the current step. Save the integrator data for this step as the solution.
                this->storage_ptr->save_data(this->t_now_ptr[0], this->y_now_ptr, this->dy_now_ptr);
            }

            if (prepare_for_next_step)
            {
                // Prep for next step
                this->t_old = this->t_now_ptr[0];
                std::memcpy(this->y_old_ptr, this->y_now_ptr, sizeof(double) * this->num_y);
                std::memcpy(this->dy_old_ptr, this->dy_now_ptr, sizeof(double) * this->num_dy);
            }
        }
    }

    // Note this is not an "else" block because the integrator may have finished with that last step.
    // Check status again to see if we are finished or there was an error in the last step
    if (this->status != 0)
    {
        // Update integration message
        this->storage_ptr->error_code = this->status;
        this->storage_ptr->success    = false;
        switch (this->status)
        {
        case 2:
            this->storage_ptr->update_message("Integration storage changed but integrator was not reset. Call `.reset()` before integrating after change.");
            break;
        case 1:
            this->storage_ptr->update_message("Integration completed without issue.");
            this->storage_ptr->success = true;
            break;
        case -1:
            this->storage_ptr->update_message("Error in step size calculation:\n\tRequired step size is less than spacing between numbers.");
            break;
        case -2:
            this->storage_ptr->update_message("Maximum number of steps (set by user) exceeded during integration.");
            break;
        case -3:
            this->storage_ptr->update_message("Maximum number of steps (set by system architecture) exceeded during integration.");
            break;
        case -4:
            this->storage_ptr->update_message("Error in step size calculation:\n\tError in step size acceptance.");
            break;
        case -9:
            this->storage_ptr->update_message("Error in CySolver initialization.");
            break;
        default:
            this->storage_ptr->update_message("Unknown status encountered during integration.");
            break;
        }
        
        // Call the finalizer on the storage class instance
        this->storage_ptr->finalize();
    }
}


void CySolverBase::change_storage(std::shared_ptr<CySolverResult> new_storage_ptr, bool auto_reset)
{   
    // Change the storage reference. 
    this->storage_ptr = new_storage_ptr;
    
    // Set status to indicate that storage has changed but reset has not been updated
    this->status = 2;

    // Changing storage requires a reset
    if (auto_reset)
    {
        this->reset();
    }
}


// Main Solve Method!
void CySolverBase::solve()
{
    while (this->check_status())
    {
        this->take_step();
    }
}

/* Dense Output Methods */
CySolverDense* CySolverBase::p_dense_output_heap()
{
    return new CySolverDense(
        this->integration_method,
        this->t_old,
        this->t_now_ptr[0],
        this->y_old_ptr,
        this->num_y,
        0 // Fake Q order just for consistent constructor call
        );
}

void CySolverBase::p_dense_output_stack(CySolverDense& dense_output_ptr)
{
    // Don't do anything. Subclasses will override this method.
}


/* PySolver Methods */
void CySolverBase::set_cython_extension_instance(PyObject* cython_extension_class_instance, DiffeqMethod py_diffeq_method)
{
    this->use_pysolver = true;
    if (cython_extension_class_instance) [[likely]]
    {
        this->cython_extension_class_instance = cython_extension_class_instance;
        this->py_diffeq_method = py_diffeq_method;

        // Change diffeq binding to the python version
        this->diffeq = &CySolverBase::py_diffeq;

        // Import the cython/python module (functionality provided by "pysolver_api.h")
        const int import_error = import_CyRK__cy__pysolver_cyhook();
        if (import_error)
        {
            this->use_pysolver = false;
            this->status = -1;
            this->storage_ptr->error_code = -51;
            this->storage_ptr->update_message("Error encountered importing python hooks.\n");
        }
        else
        {
            // TODO: Do we need to decref this at some point? During CySolver's deconstruction?
            Py_XINCREF(this->cython_extension_class_instance);
        }
    }
}

void CySolverBase::py_diffeq()
{
    // Call the differential equation in python space. Note that the optional arguments are handled by the python 
    // wrapper class. `this->args_ptr` is not used.
    call_diffeq_from_cython(this->cython_extension_class_instance, this->py_diffeq_method);
}
